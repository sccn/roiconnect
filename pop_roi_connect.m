% pop_roi_connect - call roi_connect to compute connectivity between ROIs
%
% Usage:
%  EEG = pop_roi_connect(EEG, 'key', 'val', ...);
%
% Inputs:
%  EEG - EEGLAB dataset containing ROI activity
%
% Optional inputs:
%  'trgc'      - ['on'|'off'] compute time-reverse Granger Causality. Default
%                is 'on'.
%  'crossspec' - ['on'|'off'] compute cross-spectrum from which coherence can
%                be derived. Default is 'on'.
%
% Output:
%  EEG - EEGLAB dataset with field 'roi' containing connectivity info.
%
% Note: Optional inputs to roi_connectivity_process() are also accepted.
%
% Author: Arnaud Delorme, UCSD, 2019
%
% Example
%   p = fileparts(which('eeglab')); % path
%   EEG = pop_roi_connect(EEG, 'headmodel', ...
%   EEG.dipfit.hdmfile, 'elec2mni', EEG.dipfit.coord_transform, ...
%   'sourcemodel', fullfile(p, 'functions', 'supportfiles', ...
%   'head_modelColin27_5003_Standard-10-5-Cap339.mat'), 'sourcemodel2mni', ...
%   [0 -26.6046230000 -46 0.1234625600 0 -1.5707963000 1000 1000 1000]);
%
% Use pop_roi_connectivity(EEG) to compute conectivity

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

% TO DO - Arno
% - Centralize reading head mesh and Atlas (there might be a function in
% Fieldtrip to do that) ft_read_volume ft_read_mesh
% - Make compatible with all Fieldtrip and FSL Atlases
% - Downsampling of Atlas - check bug submitted to Fieldtrip
% - Plot inside(blue) vs outside(red) voxels for source volume

function [EEG,com] = pop_roi_connect(EEG, varargin)

com = '';
if nargin < 1
    help pop_roi_connect;
    return
end

if ~isfield(EEG(1), 'roi') || ~isfield(EEG(1).roi, 'source_roi_data')
    error('Cannot find ROI data - compute ROI data first');
end

if nargin < 2
    cb_crossspec = [ 'if get(gcbo, ''value''),' ...
                     '    set(findobj(gcbf, ''userdata'', ''crossspec''), ''enable'', ''on'');' ...
                     'else,' ...
                     '    set(findobj(gcbf, ''userdata'', ''crossspec''), ''enable'', ''off'');' ...
                     'end;' ];
    cb_trgc = [ 'if get(gcbo, ''value''),' ...
                     '    set(findobj(gcbf, ''userdata'', ''trgc''), ''enable'', ''on'');' ...
                     'else,' ...
                     '    set(findobj(gcbf, ''userdata'', ''trgc''), ''enable'', ''off'');' ...
                     'end;' ];
    
    rowg = [0.1 0.6 1 0.2];
    % uigeom = { 1 1 rowg rowg 1 rowg rowg [0.1 0.6 0.9 0.3] 1 rowg 1 [0.5 1 0.35 0.5] [0.5 1 0.35 0.5] [0.5 1 0.35 0.5] [1] [0.9 1.2 1] };
    uigeom = { [1] [0.2 1 0.35] [1] [0.2 1 0.35] [0.2 1 0.35] };
    uilist = { { 'style' 'checkbox' 'string' 'Compute cross-spectrum' 'tag' 'crossspec' 'value' 1 'callback' cb_crossspec } ...
               {} { 'style' 'text' 'string' 'Bootstrap if any (n)'  'userdata' 'crossspec' }       { 'style' 'edit' 'string' '' 'tag' 'naccu1'  'userdata' 'crossspec' }    ...
               { 'style' 'checkbox' 'string' 'Compute TRGC' 'tag' 'trgc' 'value' 0 'callback' cb_trgc  } ...
               {} { 'style' 'text' 'string' 'Autoregressive model order'  'userdata' 'trgc' 'enable' 'off' } { 'style' 'edit' 'string' '20' 'tag' 'morder' 'userdata' 'trgc' 'enable' 'off'} ...
               {} { 'style' 'text' 'string' 'Bootstrap if any (n)'  'userdata' 'trgc'  'enable' 'off'}       { 'style' 'edit' 'string' '' 'tag' 'naccu2'  'userdata' 'trgc' 'enable' 'off'} };
                ...
    [result,usrdat,~,out] = inputgui('geometry', uigeom, 'uilist', uilist, 'helpcom', 'pophelp(''pop_roi_connect'')');
    if isempty(result), return, end
                      
    % check we have the same naccu
    if ~isempty(out.naccu1) && ~isempty(out.naccu2)
        if ~isequal(out.naccu1, out.naccu2)
            error('Bootstrap naccu must be the same for all methods');
        end
        if ~isempty(out.naccu1)
            out.naccu2 = out.naccu1;
        end
    end
    options = {  ...
        'morder' str2num(out.morder) ...
        'naccu' str2num(out.naccu2) ...
        'trgc'  fastif(out.trgc, 'on', 'off') ...
        'crossspec' fastif(out.crossspec, 'on', 'off') };
else
    options = varargin;
end

% process multiple datasets
% -------------------------
if length(EEG) > 1
    % check that the dipfit settings are the same
    if nargin < 2
        [ EEG, com ] = eeg_eval( 'pop_roi_connect', EEG, 'warning', 'on', 'params', options );
    else
        [ EEG, com ] = eeg_eval( 'pop_roi_connect', EEG, 'params', options );
    end
    return;
end

EEG = roi_connect(EEG, options{:});

if nargout > 1
    com = sprintf( 'EEG = pop_roi_connect(EEG, %s);', vararg2str( options ));
end

