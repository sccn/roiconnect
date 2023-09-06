% Wrapper function to compute surrogate FC metrics that are based on the
% cross-spectrum, incl. MIM/MIC, cCOH, iCOH, aCOH and wPLI.
%
% Usage: 
%   conn = shuffle_MIM(data, npcs, output, nshuf, 'freqresolution', <freqresolution>, 'roi_selection', <roi_selection>); 
%
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

function conn = shuffle_MIM(data, npcs, output, nshuf, varargin)
    % Copyright (c) 2022 Franziska Pellegrini and Stefan Haufe

    % decode input parameters
    g = finputcheck(varargin, { ...
        'freqresolution'  'integer'  { }    0;
        'roi_selection'   'cell'     { }    { }; ...
        'poolsize'        'integer'  { }     1 }, 'shuffle_MIM'); 
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

    CSpara = [];
    CSpara.subave = 0;
    CSpara.mywindow = hanning(ndat) ./ sqrt(hanning(ndat)' * hanning(ndat));
    CSpara.freqresolution = g.freqresolution;
    
%     warning('One iteration takes about 90 seconds.')
    fprintf('Generating null distribution using %d shuffles...\n', nshuf)
    fprintf('Progress of %d:\n', nshuf);

    parpool(g.poolsize)
    parfor ishuf = 1:nshuf 
        if mod(ishuf, 10) == 0
            fprintf('%d', ishuf);
        elseif mod(ishuf, 2) == 0
            fprintf('.');
        end
        
        % shuffle trials
        if ishuf == 1
            shuf_inds = 1:nepo; % true MIM values
        else
            shuf_inds = randperm(nepo);   
        end
        
        data_shuf = data(:, :, shuf_inds);
        [CS, ~, wPLI, ~] = data2cs_event_shuf(data(:, :)', data_shuf(:, :)', ndat, floor(ndat/2), ndat, [], CSpara);
%         CS = fp_tsdata_to_cpsd(data, fres, 'WELCH', 1:nchan, 1:nchan,1:nepo,shuf_inds);
        nfreqs = size(CS, 3);

        if ~isempty(intersect(output, {'MIM', 'MIC', 'cCOH' 'iCOH', 'aCOH'}))
            cCOH = zeros(nchan, nchan, nfreqs)
            aCOH = zeros(nchan, nchan, nfreqs)
            iCOH = zeros(nchan, nchan, nfreqs)
            for ifreq = 1:size(CS,3)
                pow = real(diag(CS(:,:,ifreq)));
                cCOH(:,:,ifreq) = CS(:,:,ifreq) ./ sqrt(pow*pow');
                aCOH(:,:,ifreq) = abs(CS(:,:,ifreq)) ./ sqrt(pow*pow');
                iCOH(:,:,ifreq) = abs(imag(CS(:,:,ifreq) ./ sqrt(pow*pow')));
            end
        end
        
        % loop over sender/receiver combinations to compute time-reversed GC
        MIM2 = zeros(nfreqs, ninds)
        MIC2 = zeros(nfreqs, ninds)
        for iind = 1:ninds
            if ~isequal(inds{iind}{1}, inds{iind}{2})
                %ind configuration
                subset = [inds{iind}{1} inds{iind}{2}];
                subinds = {1:length(inds{iind}{1}), length(inds{iind}{1}) + (1:length(inds{iind}{2}))};
                
                %MIC and MIM
                [MIC2(:, iind) , MIM2(:, iind)] =  roi_mim2(cCOH(subset, subset, :), subinds{1}, subinds{2});
            end
        end         
        
        % reshape (in the case of MIM) or only keep the first principal component (other metrics)
        MIM_s(:, :, :, ishuf) = get_connect_mat(MIM2, nROI, +1);
        MIC_s(:, :, :, ishuf) = get_connect_mat(MIC2, nROI, +1);
        CS_s(:, :, :, ishuf) = rm_components(permute(CS, [3 1 2 4]), npcs(1));
        cCOH_s(:, :, :, ishuf) = rm_components(permute(cCOH, [3 1 2 4]), npcs(1));
        aCOH_s(:, :, :, ishuf) = rm_components(permute(aCOH, [3 1 2 4]), npcs(1));
        iCOH_s(:, :, :, ishuf) = rm_components(permute(iCOH, [3 1 2 4]), npcs(1));
        wPLI_s(:, :, :, ishuf) = rm_components(wPLI, npcs(1));
    end

    % save methods in a struct
    clear out
    conn.inds = inds;
    for iout = 1:length(output)
        eval(['conn.' output{iout} ' = ' output{iout} '_s;'])
    end
    % shut down current parallel pool
    poolobj = gcp('nocreate');
    delete(poolobj);
    fprintf('\n');
end