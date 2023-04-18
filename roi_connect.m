% roi_connect() - compute connectivity between ROIs
%
% Usage:
%  EEG = roi_connect(EEG, 'key', 'val', ...); 
%
% Inputs:
%  EEG - EEGLAB dataset with ROI activity computed
%
% Optional inputs (choose at least one):
%  'morder'    - [integer] order of autoregressive model. Default is 20.
%  'naccu'     - [integer] number of accumulation for stats. Default is 0.
%  'methods'   - [cell] Cell of strings corresponding to methods.
%                       'CS'    : Cross spectrum
%                       'COH'   : Coherence
%                       'GC'    : Granger Causality
%                       'TRGC'  : Time-reversed Granger Causality
%                       'wPLI'  : Weighted Phase Lag Index
%                       'PDC'   : Partial directed coherence
%                       'TRPDC' : Time-reversed partial directed coherence
%                       'DTF'   : Directed transfer entropy
%                       'TRDTF' : Time-reversed directed transfer entropy
%                       'MIM'   : Multivariate Interaction Measure for each ROI
%                       'MIC'   : Maximized Imaginary Coherency for each ROI
%  'freqresolution'   - [integer] Desired frequency resolution (in number of frequencies). If
%                       specified, the signal is zero padded accordingly.
%                       Default is 0 (means no padding).
%
% Output:
%   EEG - EEG structure with EEG.roi field updated and now containing
%         connectivity information.

% Copyright (C) Arnaud Delorme, arnodelorme@gmail.com
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

function EEG = roi_connect(EEG, varargin)

    if nargin < 2
        help roi_connect;
        return
    end

    if ~isfield(EEG, 'roi') || ~isfield(EEG.roi, 'source_roi_data')
        error('Cannot find ROI data - compute ROI data first');
    else
        source_roi_data = EEG.roi.source_roi_data;
    end

    % decode input parameters
    % -----------------------
    g = finputcheck(varargin, { ...
        'morder'          'integer'  { }            20;
        'naccu'           'integer'  { }            0;
        'methods'         'cell'     { }            {}; 
        'freqresolution'  'integer'  { }            0}, 'roi_connect');    
    if ischar(g), error(g); end
    if isempty(g.naccu), g.naccu = 0; end
    tmpMethods = setdiff(g.methods, {  'CS' 'COH' 'GC' 'TRGC' 'wPLI' 'PDC' 'TRPDC' 'DTF' 'TRDTF' 'MIM' 'MIC' 'PAC'});
    if ~isempty(tmpMethods)
        error('Unknown methods %s', vararg2str(tmpMethods))
    end

    inds = {}; ninds = 0;
    nROI = EEG.roi.nROI;
    nPCA = EEG.roi.nPCA;
    for iroi = 1:nROI
        for jroi = (iroi+1):nROI
            inds{ninds+1} = {(iroi-1)*nPCA + [1:nPCA], (jroi-1)*nPCA + [1:nPCA]};
            ninds = ninds + 1;
        end
    end

    % wPLI, MIC, MIM, GC and TRGC use data2strcgmim, remaining metrics use data2spwctrgc
    methodset1 = { 'wPLI' 'MIM' 'MIC' 'GC' 'TRGC' };
    methodset2 = { 'CS' 'COH' 'PSD' 'PSDROI', 'PDC' 'TRPDC' 'DTF' 'TRDTF' };

    tmpMethods1 = intersect(g.methods, methodset1);
    if ~isempty(tmpMethods1)
        conn_mult = data2sctrgcmim(source_roi_data, EEG.pnts, g.morder, 0, g.naccu, [], inds, tmpMethods1, [], 'freqresolution', g.freqresolution);
        fields = fieldnames(conn_mult);
        for iField = 1:length(fields)
            EEG.roi.(fields{iField}) = conn_mult.(fields{iField});
        end

        % rearrange matrices
        for iMethods = 1:length(tmpMethods1)
            if strcmpi(tmpMethods1{iMethods}, 'MIM') || strcmpi(tmpMethods1{iMethods}, 'MIC')
                MI = EEG.roi.(tmpMethods1{iMethods})(:, :);
                EEG.roi.(tmpMethods1{iMethods}) = get_connect_mat( MI, EEG.roi.nROI, +1);
            elseif  strcmpi(tmpMethods1{iMethods}, 'GC') || strcmpi(tmpMethods1{iMethods}, 'TRGC')
                TRGCnet = EEG.roi.(tmpMethods1{iMethods})(:, :, 1) - EEG.roi.(tmpMethods1{iMethods})(:, :, 2);
                EEG.roi.(tmpMethods1{iMethods}) = get_connect_mat( TRGCnet, EEG.roi.nROI, -1); 
            else
                MI = EEG.roi.(tmpMethods1{iMethods})(:, :);
                EEG.roi.(tmpMethods1{iMethods}) = get_connect_mat( MI, EEG.roi.nROI, +1);
            end
        end
    end

    tmpMethods2 = intersect(g.methods, methodset2);
    if ~isempty(tmpMethods2)
        conn_mult = data2spwctrgc(source_roi_data, EEG.pnts, g.morder, 0, g.naccu, [], tmpMethods2, [], 'freqresolution', g.freqresolution);
        fields = fieldnames(conn_mult);
        for iField = 1:length(fields)
            EEG.roi.(fields{iField}) = conn_mult.(fields{iField});
        end

        % rearrange matrices
        for iMethods = 1:length(tmpMethods2)
            MI = EEG.roi.(tmpMethods2{iMethods})(:, :);
            EEG.roi.(tmpMethods2{iMethods}) = get_connect_mat( MI, EEG.roi.nROI, +1);
        end
    end

    function measure = get_connect_mat( measureOri, nROI, signVal)
        % create a ROI x ROI connectivity matrix, if needed
        % TRGCmat(f, ii, jj) is net TRGC from jj to ii
        measure = [];
        iinds = 0;
        for iroi = 1:nROI
            for jroi = (iroi+1):nROI
                iinds = iinds + 1;
                measure(:, iroi, jroi) = signVal * measureOri(:, iinds);
                measure(:, jroi, iroi) = measureOri(:, iinds);
            end
        end
    end
end