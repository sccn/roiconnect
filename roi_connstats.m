% roi_connstats() - Compute connectivity between ROIs and run statistics.
%
% Usage:
%  EEG = roi_connstats(EEG, fcomb, nshuf); 
%
% Inputs:
%  EEG           - EEGLAB dataset with ROI activity computed.
%  methods       - [cell] Cell of strings corresponding to methods.
%                       'CS'    : Cross spectrum
%                       'aCOH'  : Coherence
%                       'cCOH'  : (complex-valued) coherency
%                       'iCOH'  : absolute value of the imaginary part of coherency
%                       'GC'    : Granger Causality
%                       'TRGC'  : Time-reversed Granger Causality
%                       'wPLI'  : Weighted Phase Lag Index
%                       'PDC'   : Partial directed coherence
%                       'TRPDC' : Time-reversed partial directed coherence
%                       'DTF'   : Directed transfer entropy
%                       'TRDTF' : Time-reversed directed transfer entropy
%                       'MIM'   : Multivariate Interaction Measure for each ROI
%                       'MIC'   : Maximized Imaginary Coherency for each ROI
%  nshuf'        - [integer] number of shuffles for statistical significance testing. The first shuffle is the true value. Default is 101. 
%  roi_selection - [cell array of integers] Cell array of ROI indices {1, 2, 3, ...} indicating for which regions (ROIs) connectivity should be computed. 
%                  Default is all (set to EEG.roi.nROI).

function EEG = roi_connstats(EEG, methods, nshuf, roi_selection)

    if nargin < 2
        help roi_connstats;
        return
    end

    if ~isfield(EEG, 'roi') || ~isfield(EEG.roi, 'source_roi_data')
        error('Cannot find ROI data - compute ROI data first');
    else
        data = EEG.roi.source_roi_data;
    end

    if ~isempty(intersect(methods, {'COH'}))
        warning("'COH' is not supported anymore and will be replaced with aCOH (coherence). " + ...
            "Please double-check with the documentation if this is what you want.")
        coh_idx = strcmpi(methods, 'COH');
        g.methods{coh_idx} = 'aCOH';
    end

    if length(methods) == 1 && ~strcmpi(methods{1}, 'MIM' )
        error('Statistics are only supported for MIM, more methods will be added.')
    else
        warning('Statistics are only supported for MIM, more methods will be added.')
    end

    methodset1 = { 'MIM' }; % the remaining methods will be added eventually
    tmpMethods1 = intersect(methods, methodset1);
    if ~isempty(tmpMethods1)
        npcs = repmat(EEG.roi.nPCA, 1, EEG.roi.nROI);
        MIM_s = shuffle_MIM(data, npcs, EEG.srate, nshuf); % (nfreq, nROI, nROI, nshuf)
        EEG.roi.MIM = MIM_s;
    end
end