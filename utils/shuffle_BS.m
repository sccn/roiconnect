% Wrapper function to compute surrogate FC metrics that are based on the bispectrum, incl. PAC.
%
% Usage: 
%   conn = shuffle_BS(data, npcs, output, nshuf, fs); (without optional inputs) 
%   conn = shuffle_BS(data, npcs, output, nshuf, fs, 'freqresolution', <freqresolution>, 'roi_selection', <roi_selection>, 'poolsize', <poolsize>); 
%
% Inputs:
%   data             - (nchan x len_epoch x ntrials) source ROI data
%   npcs             - (1 x nROIs) vector, each entry contains the number of principal components (PCs)
%   output           - [cell array of string] Cell array of methods e.g. {'CS' 'MIM' 'wPLI' 'cCOH' 'aCOH' 'iCOH'}
%   nshuf            - [integer] number of shuffles
%   fs               - [integer] sampling rate in Hz
%   fcomb            - [struct] Frequency combination for which PAC is computed (in Hz). Must have fields 'low' and 
%                      'high' with fcomb.low < fcomb.high. For example, fcomb.low = 10 and fcomb.high = 50 if single 
%                      frequencies are used. fcomb.low = [4 8] and fcomb.high = [48 50] if frequency bands are used 
%                      (might take a long time to compute, so use with caution). Default is {} (this will cause an error).
% 
% Optional inputs:
%   'freqresolution' - [integer] Desired frequency resolution (in number of frequencies). If
%                       specified, the signal is zero padded accordingly. Default is 0 (means no padding).
%   'roi_selection'  - [cell array of integers] Cell array of ROI indices {1, 2, 3, ...} indicating for which regions (ROIs) connectivity should be computed. 
%                       Default is all (set to EEG.roi.nROI).
%   'poolsize'       - [integer] Number of workers in the parallel pool (check parpool documentation) for parallel computing
% 
% Outputs
%   conn             - [struct] Struct of (nfreq x nROI x nROI x nshuf) FC metrics
%
% Authors: 
%   Zixuan Liu, zixuan.liu@campus.tu-berlin.de
%   Tien Dung Nguyen, tien-dung.nguyen@charite.de 

function conn = shuffle_BS(data, npcs, output, nshuf, fs, fcomb, varargin)

    % decode input parameters
    g = finputcheck(varargin, { ...
        'freqresolution'  'integer'  { }    0;
        'roi_selection'   'cell'     { }    { }; ...
        'poolsize'        'integer'  { }     1 ;
        }, 'shuffle_BS'); 
    if ischar(g), error(g); end

    if ~isfield(fcomb, 'low') || ~isfield(fcomb, 'high')
        help roi_pac;
        error('Frequency pair cannot be found - check the documentation for the "fcomb" input parameter in shuffle_BS.')
    end

    if fcomb.high < fcomb.low
        help roi_pac;
        error('fcomb.high must be smaller than fcomb.low - check the documentation for the "fcomb" input parameter in "shuffle_BS".')
    end
    
    [nchan, ndat, ~] = size(data);

    % [inds, PCA_inds] = fp_npcs2inds(npcs);
    if isempty(g.roi_selection)
        nROI = nchan/npcs(1);
    else
        nROI = length(g.roi_selection);
    end
    nPCA = npcs(1);

    % choose ROIs if desired, take number of PCs into account
    if ~isempty(g.roi_selection)
        data_new = zeros(nROI * nPCA, size(data, 2), size(data, 3));
        
        start_idx_new  = 1;
        end_idx_new = nPCA;
        for iroi = 1:nROI
            end_idx = g.roi_selection{iroi} * nPCA;
            start_idx = end_idx - (nPCA - 1);
            data_new(start_idx_new:end_idx_new, :, :) = data(start_idx:end_idx, :, :);
    
            start_idx_new = start_idx_new + nPCA;
            end_idx_new = start_idx_new + nPCA - 1;
        end
        data = data_new;
    end

    % only keeep first PC
    if nPCA > 1
        warning('Only the first principal component will be used to determine PAC.')
        data = data(1:nPCA:end, :, :);
    end
    
    % warning('One iteration takes about 90 seconds.')
    fprintf('Generating null distribution using %d shuffles...\n', nshuf)
    fprintf('Progress of %d:\n', nshuf);
    
    % Define bispectrum parameters
    fres = fs;
    frqs = sfreqs(fres, fs);
    
    % extract all individual frequencies in the selected bands
    size_low = size(fcomb.low, 2);
    size_high = size(fcomb.high, 2);
    mask_inds_low = frqs >= fcomb.low(1) & frqs <= fcomb.low(size_low);
    mask_inds_high = frqs >= fcomb.high(1) & frqs <= fcomb.high(size_high);
    frqs_low = frqs(mask_inds_low); 
    frqs_high = frqs(mask_inds_high);
    
    % determine all frequency combinations
    [m, n] = ndgrid(frqs_low, frqs_high);
    frqs_combs = [m(:),n(:)]; 
    n_combs = size(frqs_combs, 1);
    if n_combs > 20
        % according to our test simulations, the computation time scales linearly with the number of frequency pairs times 2, assuming no other ongoing CPU-heavy processes
        time_est = 2 * n_combs; 
        warning('PAC is going to be estimated on %d frequency pair(s). Estimated time: %d seconds', n_combs, time_est);
    end
    
    freqinds_low = zeros(n_combs, 2);
    freqinds_up = zeros(n_combs, 2);
    for i = 1:n_combs
        low = frqs_combs(i,1);
        high = frqs_combs(i,2);
        freqinds_low(i,:) = [find(frqs == low) find(frqs == high - low)]; 
        freqinds_up(i,:) = [find(frqs == low) find(frqs == high)];
    end

    % check if Parallel Processing Toolbox is available and licensed
    if license('test', 'Distrib_Computing_Toolbox') && ~isempty(ver('parallel'))
        if isfield(g, 'poolsize') && isnumeric(g.poolsize) && g.poolsize > 0
            % check if there's already an existing parallel pool
            currentPool = gcp('nocreate');
            if isempty(currentPool)
                parpool(g.poolsize);
            end
        end
    else
        disp('Parallel Processing Toolbox is not installed or licensed.');
    end

    % Pre-allocate variables outside the parfor loop
    b_orig = zeros(nROI, nROI, nshuf);
    b_anti = zeros(nROI, nROI, nshuf);  
    % b_orig_norm = zeros(nROI, nROI, nshuf);
    % b_anti_norm = zeros(nROI, nROI, nshuf);

    parfor ishuf = 1:nshuf
        b_orig_ishuf = zeros(nROI, nROI);
        b_anti_ishuf = zeros(nROI, nROI);
        % b_orig_norm_ishuf = zeros(nROI, nROI);
        % b_anti_norm_ishuf = zeros(nROI, nROI);

        % Iterate over ROI pairs
        for proi = 1:nROI
            for aroi = proi:nROI
                % Compute bispectrum % nchan by nchan by nchan by number_of_peaks by number_of_shuffles 
                X = data([proi aroi],:,:); % number of regions x epoch length x ishuf
                [BS, ~] = fp_data2bs_event_uni(X(:, :)', ndat, floor(ndat/2), ndat, freqinds_up, ishuf); % pass (f1,f2) through freqinds_up
                BS_low = BS(:,:,:,1);
                BS_up = BS(:,:,:,2);
    
                [RTP_low,~] = data2bs_threenorm(X(:, :)', ndat, floor(ndat/2), ndat, freqinds_low);
                [RTP_up,~] = data2bs_threenorm(X(:, :)', ndat, floor(ndat/2), ndat, freqinds_up);
    
                [biv_orig_low, biv_anti_low] = calc_pac(BS_low, RTP_low); 
                [biv_orig_up, biv_anti_up] = calc_pac(BS_up, RTP_up);
    
                % PAC_km(f1, f2) = 0.5 * |Bkmm(f1, f2-f1)| + 0.5 * |Bkmm(f1, f2)|
                b_orig_ishuf(aroi,proi) = mean([biv_orig_up(1) biv_orig_low(1)]); 
                b_orig_ishuf(proi,aroi) = mean([biv_orig_up(2) biv_orig_low(2)]);
                b_anti_ishuf(aroi,proi) = mean([biv_anti_up(1) biv_anti_low(1)]);  
                b_anti_ishuf(proi,aroi) = mean([biv_anti_up(2) biv_anti_low(2)]); 
    
                % % normalized versions (for bicoherence)
                % b_orig_norm(aroi,proi) = mean([biv_orig_up_norm(1) biv_orig_low_norm(1)]);
                % b_orig_norm(proi,aroi) = mean([biv_orig_up_norm(2) biv_orig_low_norm(2)]);
                % b_anti_norm(aroi,proi) = mean([biv_anti_up_norm(1) biv_anti_low_norm(1)]);  
                % b_anti_norm(proi,aroi) = mean([biv_anti_up_norm(2) biv_anti_low_norm(2)]);
            end
            
        end

        % Store the shuffle information
        b_orig(:,:, ishuf) = b_orig_ishuf;
        b_anti(:,:, ishuf) = b_anti_ishuf;
        % b_orig_norm(:,:, ishuf) = b_orig_norm_ishuf;
        % b_anti_norm(:,:, ishuf) = b_anti_norm_ishuf;

    end


    clear out

    % Save PAC results in the output structure
    conn.PAC.b_orig = b_orig;
    conn.PAC.b_anti = b_anti;
    % conn.PAC.b_orig_norm = b_orig_norm;
    % conn.PAC.b_anti_norm = b_anti_norm;

    % shut down current parallel pool only if the toolbox is available
    if license('test', 'Distrib_Computing_Toolbox') && ~isempty(ver('parallel'))
        poolobj = gcp('nocreate');
        if ~isempty(poolobj)
            delete(poolobj);
        end
    end
    fprintf('\n');
end