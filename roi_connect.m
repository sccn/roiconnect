% roi_connect() - compute connectivity between ROIs
%
% Usage:
%  EEG = roi_connect(EEG, 'key', 'val', ...); 
%
% Inputs:
%  EEG - EEGLAB dataset with ROI activity computed
%
% Optional inputs (choose at least one):
%  'trgc'      - ['on'|'off'] compute time-reverse Granger Causality. Default
%                is 'on'.
%  'crossspec' - ['on'|'off'] compute cross-spectrum from which coherence can
%                be derived. Default is 'on'.
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
    'morder'      'integer' { }            20;
    'naccu'       'integer' { }            0;
    'methods'     'cell'    { }            {} }, 'roi_connect');    
if ischar(g), error(g); end
if isempty(g.naccu), g.naccu = 0; end

inds = {}; ninds = 0;
nROI = EEG.roi.nROI;
nPCA = EEG.roi.nPCA;
for iroi = 1:nROI
    for jroi = (iroi+1):nROI
        inds{ninds+1} = {(iroi-1)*nPCA + [1:nPCA], (jroi-1)*nPCA + [1:nPCA]};
        inds{ninds+2} = {(jroi-1)*nPCA + [1:nPCA], (iroi-1)*nPCA + [1:nPCA]};
        ninds = ninds + 2;
    end
end
if ~isempty(g.methods)
    conn_mult = data2sctrgcmim(source_roi_data, EEG.srate, g.morder, 0, g.naccu, [], inds, g.methods);
    
    % can compute CS below as well (not sure if it is better than the one
    % above, but it is different
    % conn_uni  = data2spwctrgc2(source_roi_data, EEG.srate, g.morder, 0, g.naccu, [], {'CS'});
end

fields = fieldnames(conn_mult);
for iField = 1:length(fields)
    EEG.roi.(fields{iField}) = conn_mult.(fields{iField});
end
