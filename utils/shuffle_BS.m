% Wrapper function to compute surrogate FC metrics that are based on the bispectrum, incl. PAC.
%
% Usage: 
%   conn = shuffle_BS(data, npcs, output, nshuf, 'freqresolution', <freqresolution>, 'roi_selection', <roi_selection>); 
%...............
% Inputs:
%   data             - (nchan x len_epoch x ntrials) source ROI data
%   npcs             - (1 x nROIs) vector, each entry contains the number of principal components (PCs)
%   output           - [cell array of string] Cell array of methods e.g. {'CS' 'MIM' 'wPLI' 'cCOH' 'aCOH' 'iCOH'}
%   nshuf            - [integer] number of shuffles
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

function conn = shuffle_BS(data, npcs, output, nshuf, fs, varargin)
    % Copyright (c) 2022 Franziska Pellegrini and Stefan Haufe

    % decode input parameters
    g = finputcheck(varargin, { ...
        'freqresolution'  'integer'  { }    0;
        'roi_selection'   'cell'     { }    { }; ...
        'poolsize'        'integer'  { }     1 ;...
        'fcomb'           'struct'   { }    struct;
        }, 'shuffle_BS'); 
    if ischar(g), error(g); end
    
    [nchan, ndat, nepo] = size(data);

%     [inds, PCA_inds] = fp_npcs2inds(npcs);
%     ninds = length(inds);
    inds = {}; ninds = 0;
    if isempty(g.roi_selection)
        nROI = nchan/npcs(1);
    else
        nROI = length(g.roi_selection);
    end
    nPCA = npcs(1);
    for iroi = 1:nROI
        for jroi = (iroi+1):nROI
            inds{ninds+1} = {(iroi-1)*nPCA + [1:nPCA], (jroi-1)*nPCA + [1:nPCA]};
            ninds = ninds + 1;
        end
    end

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
    
    
    % Check if Parallel Processing Toolbox is available and licensed
    if license('test', 'Distrib_Computing_Toolbox') && ~isempty(ver('parallel'))
        % If g.poolsize is defined, create a parallel pool of that size
        if isfield(g, 'poolsize') && isnumeric(g.poolsize) && g.poolsize > 0
            % Check if there's already an existing parallel pool
            currentPool = gcp('nocreate');
            if isempty(currentPool)
                parpool(g.poolsize);
            else
                % Optionally, adjust current pool size to g.poolsize
                % If this is needed, delete the currentPool and then initiate a new one
                % delete(currentPool);
                % parpool(g.poolsize);
            end
        end
    else
        disp('Parallel Processing Toolbox is not installed or licensed.');
    end
    

    % Define bispectrum parameters
    fcomb = g.fcomb;
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

    % Initialize variables to store the PAC results
    PAC_orig = zeros(nROI, nROI, nshuf);
    PAC_anti = zeros(nROI, nROI, nshuf);
    PAC_orig_norm = zeros(nROI, nROI, nshuf);
    PAC_anti_norm = zeros(nROI, nROI, nshuf);

    % Iterate over ROI pairs
    for proi = 1:nROI
        for aroi = proi:nROI
            % Compute bispectrum % nchan by nchan by nchan by number_of_peaks by number_of_shuffles 
            X = data([proi aroi],:,:); % number of regions X epoch length X trails
            [BS, ~] = fp_data2bs_event_uni(X(:, :)', ndat, floor(ndat/2), ndat, freqinds_up, nshuf); % pass (f1,f2) through freqinds_up
            BS_low = BS(:,:,:,1,:);
            BS_up = BS(:,:,:,2,:);
            % Call bs2pac function
            [RTP_low,~] = data2bs_threenorm(X(:, :)', ndat, floor(ndat/2), ndat, freqinds_low);
            [RTP_up,~] = data2bs_threenorm(X(:, :)', ndat, floor(ndat/2), ndat, freqinds_up);
            
            [biv_orig_low, biv_anti_low, biv_orig_low_norm, biv_anti_low_norm] = calc_pac(BS_low, RTP_low); % add dimension
            [biv_orig_up, biv_anti_up, biv_orig_up_norm, biv_anti_up_norm] = calc_pac(BS_up, RTP_up);
            
            % PAC_km(f1, f2) = 0.5 * |Bkmm(f1, f2-f1)| + 0.5 * |Bkmm(f1, f2)|
            b_orig(aroi,proi) = mean([biv_orig_up(1) biv_orig_low(1)]); 
            b_orig(proi,aroi) = mean([biv_orig_up(2) biv_orig_low(2)]);
            b_anti(aroi,proi) = mean([biv_anti_up(1) biv_anti_low(1)]);  
            b_anti(proi,aroi) = mean([biv_anti_up(2) biv_anti_low(2)]); 
            
            % normalized versions
            b_orig_norm(aroi,proi) = mean([biv_orig_up_norm(1) biv_orig_low_norm(1)]);
            b_orig_norm(proi,aroi) = mean([biv_orig_up_norm(2) biv_orig_low_norm(2)]);
            b_anti_norm(aroi,proi) = mean([biv_anti_up_norm(1) biv_anti_low_norm(1)]);  
            b_anti_norm(proi,aroi) = mean([biv_anti_up_norm(2) biv_anti_low_norm(2)]);
    
            % Store PAC results
            PAC_orig = b_orig;
            PAC_anti = b_anti;
            PAC_orig_norm = b_orig_norm;
            PAC_anti_norm = b_anti_norm;
        end
    end

    clear out

    % Save PAC results in the output structure
    conn.PAC.b_orig = PAC_orig;
    conn.PAC.b_anti = PAC_anti;
    conn.PAC.b_orig_norm = PAC_orig_norm;
    conn.PAC.b_anti_norm = PAC_anti_norm;

    % shut down current parallel pool only if the toolbox is available
    if license('test', 'Distrib_Computing_Toolbox') && ~isempty(ver('parallel'))
        poolobj = gcp('nocreate');
        if ~isempty(poolobj)
            delete(poolobj);
        end
    end
    fprintf('\n');
end