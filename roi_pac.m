% roi_pac() - Compute phase-amplitude coupling (PAC) between ROIs  (cf. Zandvoort and Nolte, 2021). 
% Wrapper function for pac_bispec().
%
% Usage:
%  EEG = roi_pac(EEG, fcomb, bs_outopts); 
%
% Inputs:
%  EEG        - EEGLAB dataset with ROI activity computed.
%  fcomb      - [struct] Frequency combination for which PAC is computed. Must have fields
%               'low' and 'high' with fcomb.low < fcomb.high. For example, fcomb.low =
%               10 (Hz), fcomb.high = 50 (Hz).
%  bs_outopts - [integer] Option which bispectral tensors should be stored in EEG.roi.PAC. Default is 1.
%                          1 - store all tensors: b_orig, b_anti, b_orig_norm, b_anti_norm
%                          2 - only store: b_orig, b_anti
%                          3 - only store: b_orig_norm, b_anti_norm
%                          4 - only store: b_orig, b_orig_norm
%                          5 - only store: b_anti, b_anti_norm
%
% Output:
%   EEG - EEG structure with EEG.roi field updated and now containing
%         connectivity information.

function EEG = roi_pac(EEG, fcomb, bs_outopts)

    if nargin < 2
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
        error('Frequency pair cannot be found - check the documentation for the fcomb input parameter')
    end

    if fcomb.high < fcomb.low
        help roi_pac;
        error('fcomb.high must be smaller than fcomb.low - check the documentation for the fcomb input parameter')
    end

    [~, ndat, ~] = size(data);
    nROI = EEG.roi.nROI;
    segleng = ndat;
    segshift = floor(ndat/2);
    epleng = ndat;
    fs = EEG.srate;

    params = [];
    params.nROI = nROI;
    params.segleng = segleng;
    params.segshift = segshift;
    params.epleng = epleng;
    params.fcomb = fcomb;
    params.fs = fs;
    
    [b_orig, b_anti, b_orig_norm,b_anti_norm] = pac_bispec(data, params);
    
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

