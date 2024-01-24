% roi_connstats() - Compute connectivity between ROIs and run statistics.
%
% Usage:
%  EEG = roi_connstats(EEG, fcomb, nshuf); 
%
% Inputs:
%  EEG                - EEGLAB dataset with ROI activity computed.
%
% Optional inputs:
%  'methods'          - [cell] Cell of strings corresponding to methods.
%                           'CS'    : Cross spectrum
%                           'aCOH'  : Coherence
%                           'cCOH'  : (complex-valued) coherency
%                           'iCOH'  : absolute value of the imaginary part of coherency
%                           'wPLI'  : Weighted Phase Lag Index
%                           'MIM'   : Multivariate Interaction Measure for each ROI
%                           'MIC'   : Maximized Imaginary Coherency for each ROI
%  'nshuf'            - [integer] Number of shuffles for statistical significance testing. The first shuffle is the true value. Default is 1001. 
%  'roi_selection'    - [cell array of integers] Cell array of ROI indices {1, 2, 3, ...} indicating for which regions (ROIs) connectivity should be computed. 
%                       Default is all (set to EEG.roi.nROI).
%  'freqresolution'   - [integer] Desired frequency resolution (in number of frequencies). If
%                       specified, the signal is zero padded accordingly.
%                       Default is 0 (means no padding).
%   'poolsize'        - [integer] Number of workers in the parallel pool (check parpool documentation) for parallel computing
% Authors: 
%   Tien Dung Nguyen, tien-dung.nguyen@charite.de
%   Zixuan Liu, zixuan.liu@campus.tu-berlin.de

function EEG = roi_connstats(EEG, varargin)

    if nargin < 2
        help roi_connstats;
        % TODO: open GUI as in roi_connect.m
        return
    end

    if ~isfield(EEG, 'roi') || ~isfield(EEG.roi, 'source_roi_data')
        error('Cannot find ROI data - compute ROI data first');
    else
        data = EEG.roi.source_roi_data;
    end

    % decode input parameters
    % Add fcomb in parameter list
    g = finputcheck(varargin, { ...
        'methods'         'cell'     { }            { }; 
        'nshuf'           'integer'  { }            1001;
        'freqresolution'  'integer'  { }            0;
        'roi_selection'   'cell'     { }            { }; ...
        'poolsize'        'integer'  { }              1 ; ...
        'fcomb'           'struct'   { }            struct; ...
        },'roi_connstats');    
    if ischar(g), error(g); end

    if ~isempty(intersect(g.methods, {'COH'}))
        warning("'COH' is not supported anymore and will be replaced with aCOH (coherence). " + ...
            "Please double-check with the documentation if this is what you want.")
        coh_idx = strcmpi(g.methods, 'COH');
        g.methods{coh_idx} = 'aCOH';
    end

    methodset1 = { 'CS' 'MIM' 'MIC' 'wPLI' 'cCOH' 'aCOH' 'iCOH' }; % GC/TRGC, PDC/TRPDC, DTF/TRDTF not included (yet)
    methodset2 = {'PAC'} ;
    tmpMethods1 = intersect(g.methods, methodset1);
    tmpMethods2 = intersect(g.methods, methodset2);

    if ~isempty(tmpMethods1)
        npcs = repmat(EEG.roi.nPCA, 1, EEG.roi.nROI);
        conn = shuffle_CS(data, npcs, tmpMethods1, g.nshuf, 'freqresolution', g.freqresolution, 'roi_selection', g.roi_selection, 'poolsize', g.poolsize); % (nfreq, nROI, nROI, nshuf)
        for iMethod = 1:length(tmpMethods1)
            EEG.roi.(tmpMethods1{iMethod}) = conn.(tmpMethods1{iMethod});
            if strcmpi(tmpMethods1{iMethod}, 'MIM') || strcmpi(tmpMethods1{iMethod}, 'MIC')
                EEG.roi.inds = conn.inds;
            end
        end
    end
    %
    if ~isempty(tmpMethods2)
        npcs = repmat(EEG.roi.nPCA, 1, EEG.roi.nROI);
        
        fs = EEG.roi.srate;
        conn = shuffle_BS(data, npcs, tmpMethods2, g.nshuf, fs, 'freqresolution', g.freqresolution, 'roi_selection', g.roi_selection, 'poolsize', g.poolsize, 'fcomb', g.fcomb);
        for iMethod = 1:length(tmpMethods2)
            EEG.roi.(tmpMethods2{iMethod}) = conn.(tmpMethods2{iMethod});
        end
    end
end