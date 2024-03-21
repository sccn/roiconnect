% roi_pac() - Computes phase-amplitude coupling (PAC) between ROIs  (cf. Zandvoort and Nolte, 2021). 
% Wrapper function for pac_bispec().
%
% Usage:
%  EEG = roi_pac(EEG, fcomb, bs_outopts); 
%
% Inputs:
%   EEG           - EEGLAB dataset with ROI activity computed.
%   fcomb         - [struct] Frequency combination for which PAC is computed (in Hz). Must have fields 'low' and 
%                   'high' with fcomb.low < fcomb.high. For example, fcomb.low = 10 and fcomb.high = 50 if single 
%                   frequencies are used. fcomb.low = [4 8] and fcomb.high = [48 50] if frequency bands are used 
%                   (might take a long time to compute, so use with caution). Default is {} (this will cause an error).
%   bs_outopts    - [integer] Option which bispectral tensors should be stored in EEG.roi.PAC. Default is 1. "orig" means 
%                   original bispectrum, "anti" means antisymmetrized bispectrum, "norm" means normalized bispectrum.
%                          1 - store all tensors: b_orig, b_anti, b_orig_norm, b_anti_norm
%                          2 - only store: b_orig, b_anti
%                          3 - only store: b_orig_norm, b_anti_norm
%                          4 - only store: b_orig, b_orig_norm
%                          5 - only store: b_anti, b_anti_norm
%   roi_selection - [cell array of integers] Cell array of ROI indices {1, 2, 3, ...} indicating for which regions (ROIs) connectivity should be computed. 
%                  Default is all (set to EEG.roi.nROI).
%
% Output:
%   EEG - EEG structure with EEG.roi field updated and now containing
%         connectivity information.

function EEG = roi_pac(EEG, fcomb, bs_outopts, roi_selection)

    if nargin < 4
        help roi_pac;
        return
    end

    if ~isfield(EEG, 'roi') || ~isfield(EEG.roi, 'source_roi_data')
        error('Cannot find ROI data - compute ROI data first');
    else
        data = EEG.roi.source_roi_data;
    end

    if ~isfield(fcomb, 'low') || ~isfield(fcomb, 'high')
        help roi_pac;
        error('Frequency pair cannot be found - check the documentation for the "fcomb" input parameter in "roi_pac.m".')
    end

    if fcomb.high < fcomb.low
        help roi_pac;
        error('fcomb.high must be smaller than fcomb.low - check the documentation for the fcomb input parameter in "roi_pac".')
    end

    [~, ndat, ~] = size(data);
    nROI = EEG.roi.nROI;
    nPCA = EEG.roi.nPCA;
    segleng = ndat;
    segshift = floor(ndat/2);
    epleng = ndat;
    fs = EEG.roi.srate;

    params = [];
    params.nROI = nROI;
    params.segleng = segleng;
    params.segshift = segshift;
    params.epleng = epleng;
    params.fcomb = fcomb;
    params.fs = fs;
    params.roi_selection = roi_selection;
    
    % only keeep first PC
    if nPCA > 1
        warning('Only the first principal component will be used to determine PAC.')
        data = data(1:nPCA:end, :, :);
    end

    % choose ROIs if desired
    if ~isempty(roi_selection)
        nROI = length(roi_selection);
        data_new = zeros(nROI, size(data, 2), size(data, 3));
        for iroi = 1:nROI
            data_new(iroi, :, :) = data(roi_selection{iroi}, :, :);
        end
        data = data_new;
        params.nROI = nROI;
        EEG.roi.PAC.roi_selection = roi_selection;
    end
    [b_orig, b_anti, b_orig_norm, b_anti_norm] = data2bs_pac(data, params); 
    
    % options which bispectral tensors to store
    switch bs_outopts
        case 2
            EEG.roi.PAC.b_orig = b_orig;
            EEG.roi.PAC.b_anti = b_anti;
        case 3
            EEG.roi.PAC.b_orig_norm = b_orig_norm;
            EEG.roi.PAC.b_anti_norm = b_anti_norm;
        case 4
            EEG.roi.PAC.b_orig = b_orig;
            EEG.roi.PAC.b_orig_norm = b_orig_norm;
        case 5
            EEG.roi.PAC.b_anti = b_anti;
            EEG.roi.PAC.b_anti_norm = b_anti_norm;
        otherwise
            EEG.roi.PAC.b_orig = b_orig;
            EEG.roi.PAC.b_anti = b_anti;
            EEG.roi.PAC.b_orig_norm = b_orig_norm;
            EEG.roi.PAC.b_anti_norm = b_anti_norm;
    end
end

