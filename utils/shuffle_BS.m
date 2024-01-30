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
%   Franziska Pellegrini, franziska.pellegrini@charite.de
%   Stefan Haufe, haufe@tu-berlin.de
%   Tien Dung Nguyen, tien-dung.nguyen@charite.de

function conn = shuffle_BS(data, npcs, output, nshuf, varargin)
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

    % CSpara = [];
    % CSpara.subave = 0;
    % CSpara.mywindow = hanning(ndat) ./ sqrt(hanning(ndat)' * hanning(ndat));
    % CSpara.freqresolution = g.freqresolution;
    % CSpara.nshuf = g.nshuf;
    
    % warning('One iteration takes about 90 seconds.')
    fprintf('Generating null distribution using %d shuffles...\n', nshuf)
    fprintf('Progress of %d:\n', nshuf);

    % if Parallel Processing Toolbox is installed then do line 83-84 else
    % do line 85 as for not parfor
    % if Parallel Processing Toolbox is installed do parpool else not
    
    %parpool(g.poolsize)
    
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
    
    % Don;t need for loop in shuffle_BS
    % parfor ishuf = 1:nshuf 
    %     if mod(ishuf, 10) == 0
    %         fprintf('%d', ishuf);
    %     elseif mod(ishuf, 2) == 0
    %         fprintf('.');
    %     end
        
        % shuffle trials
        % if ishuf == 1
        %     shuf_inds = 1:nepo; % true MIM values
        % else
        %     shuf_inds = randperm(nepo);   
        % end
        % 
        % data_shuf = data(:, :, shuf_inds);
    
        % Starts here
        %[CS, ~, wPLI, ~] = data2cs_event_shuf(data(:, :)', data_shuf(:, :)', ndat, floor(ndat/2), ndat, [], CSpara); 
        % CS = fp_tsdata_to_cpsd(data, fres, 'WELCH', 1:nchan, 1:nchan,1:nepo,shuf_inds);

    % Define bispectrum parameters
    fcomb = g.fcomb;
    if isstruct(fcomb)
        fcomb = [fcomb.low, fcomb.high];
    else
        fcomb = fcomb;
    end
    
    % Compute bispectrum % nchan by nchan by nchan by number_of_peaks by number_of_shuffles 
    [BS, RTP] = fp_data2bs_event_uni(data(:, :)', ndat, floor(ndat/2), ndat, fcomb, nshuf); 
    
    % Initialize variables to store the PAC results
    PAC_orig = zeros(nROI, nROI, nshuf);
    PAC_anti = zeros(nROI, nROI, nshuf);
    PAC_orig_norm = zeros(nROI, nROI, nshuf);
    PAC_anti_norm = zeros(nROI, nROI, nshuf);
    
    % Iterate over ROI pairs
    for proi = 1:nROI
        for aroi = proi:nROI
            % Call bs2pac function
            [biv_orig, biv_anti, biv_orig_norm, biv_anti_norm] = calc_pac(BS(proi, aroi, :, :), RTP(proi, aroi, :, :));
    
            % Store PAC results
            PAC_orig(proi, aroi, :) = biv_orig;
            PAC_anti(proi, aroi, :) = biv_anti;
            PAC_orig_norm(proi, aroi, :) = biv_orig_norm;
            PAC_anti_norm(proi, aroi, :) = biv_anti_norm;
        end
    end

    % PAC = zeros(nfreqs, ninds);
    
    % PAC_s(:, :, :, nshuf) = get_connect_mat(PAC, nROI, +1);
    % 
    % %MIM_s(:, :, :, ishuf) = get_connect_mat(MIM2, nROI, +1);
    % %CS_s(:, :, :, ishuf) = rm_components(permute(CS, [3 1 2 4]), npcs(1));
    % %disp(CS);
    % %disp(nfreqs)
    % save methods in a struct
    clear out
    % conn.inds = inds;
    % disp(conn)
    % for iout = 1:length(output)
    %     eval(['conn.' output{iout} ' = ' output{iout} '_s;'])
    % end

    % Save PAC results in the output structure
    conn.PAC_orig = PAC_orig;
    conn.PAC_anti = PAC_anti;
    conn.PAC_orig_norm = PAC_orig_norm;
    conn.PAC_anti_norm = PAC_anti_norm;


    % shut down current parallel pool only if the toolbox is available
    if license('test', 'Distrib_Computing_Toolbox') && ~isempty(ver('parallel'))
        poolobj = gcp('nocreate');
        if ~isempty(poolobj)
            delete(poolobj);
        end
    end
    fprintf('\n');
end