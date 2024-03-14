% POP_ROI_ACTIVITY - call roi_activity to compute activities of ROIs using
%                    eLoreta
%
% Usage:
%  EEG = pop_roi_activity(EEG, 'key', 'val', ...);
%
% Inputs:
%  EEG - EEGLAB dataset. If continuous data epochs are extracted.
%
% Required inputs:
%  'headmodel'   - [string] head model file in MNI space
%  'sourcemodel' - [string] source model file
%
% Optional inputs:
%  'epochlen'         - [float] epoch length (default is 2 seconds). ROIconnect
%                       has not been tested with other epoch lenghts.
%  'epochrecur'       - [float] epoch reccurence (in second). See
%                       eeg_regepochs for more information.
%  'resample'         - [integer] resample to the desired sampling rate. Default
%                       is 100. Adujst the model order accordingly. ROIconnect
%                       has only be tested with 100 Hz sampling rates.
%  'effectchanges'    - ['on'|'off'] apply resampling and epoching to the
%                       dataset given as input. Otherwise, only update
%                       EEG.roi structure. Default is 'off'.
%  'fooof'            - ['on'|'off'] enable FOOOF analysis (this method can be 
%                       used to parameterize neural power spectra and is described here: 
%                       https://fooof-tools.github.io/fooof/). Default is 'off'.
%  'fooof_frange'     - [ ] FOOOF fitting range. Default is [1 30] like in the MATLAB example: 
%                       https://github.com/fooof-tools/fooof_mat/blob/main/examples/fooof_example_one_spectrum.m.
%  'freqresolution'   - [integer] Desired frequency resolution (in number of frequencies). If
%                       specified, the signal is zero padded accordingly.
%                       Default is 0 (means no padding).
%  'chansel'          - [cell array of string] channel selection. Default is all.
%  'lowmemory'        - ['on'|'off'] Option to run the code with low memory, though, it might take significantly longer to complete. When turned on, the estimation of voxel-wise spectral power 
%                       will require less memory.
%
% Other optional inputs:
%  All ROI_ACTIVITY parameters are accepted as input and passed on.
%
% Output:
%  EEG - EEGLAB dataset with field 'roi' containing connectivity info.
%
% Example
%   p = fileparts(which('eeglab')); % path
%   EEG = pop_roi_activity(EEG, 'headmodel', ...
%   EEG.dipfit.hdmfile, 'elec2mni', EEG.dipfit.coord_transform, ...
%   'sourcemodel', fullfile(p, 'functions', 'supportfiles', ...
%   'head_modelColin27_5003_Standard-10-5-Cap339.mat'), 'sourcemodel2mni', ...
%   [0 -26.6046230000 -46 0.1234625600 0 -1.5707963000 1000 1000 1000]);
%
% Author: Arnaud Delorme, UCSD, 2019
%
% See also ROI_ACTIVITY

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

function [EEG,com] = pop_roi_activity(EEG, varargin)

com = '';
if nargin < 1
    help pop_roi_activity;
    return
end

if ~isstruct(EEG)
    
    fig = varargin{1};
    userdata = get(fig, 'userdata');
    dipfitFile = userdata{1};
    objtmp = findobj(fig, 'tag', 'leadfieldselect');
    str1 = get(objtmp, 'string');
    str1 = str1{1};
    
    switch EEG
        case 'cb_select'
            if ~isempty(strfind(str1, 'none'))
                set(objtmp, 'value', 2);
            end
            if get(objtmp, 'value') == 1
                set(findobj(fig, 'tag', 'leadfield'), 'string', dipfitFile, 'enable', 'off');
                set(findobj(fig, 'tag', 'but'), 'enable', 'off');
            else
                set(findobj(fig, 'tag', 'leadfield'), 'string', '', 'enable', 'off');
                set(findobj(fig, 'tag', 'leadfield'), 'enable', 'on');
                set(findobj(fig, 'tag', 'but'), 'enable', 'on');
            end
        case 'cb_load'
            [tmpfilename, tmpfilepath] = uigetfile('*', 'Select a text file');
            if tmpfilename(1) ~=0, set(findobj('parent', fig, 'tag', 'leadfield'), 'string', fullfile(tmpfilepath,tmpfilename)); end
    end
    
    % update Atlas list
    leadFieldFile = get(findobj(fig, 'tag', 'leadfield'), 'string');
    iVal = 1;
    if ~isempty(leadFieldFile)
        try
            tmp = load('-mat', leadFieldFile);
            if isfield(tmp, 'atlas')
                atlasList = { 'Desikan-Kilianny' };
            else
                atlasList = { tmp.Atlas.Name };
            end
        catch
            disp('Error reading Atlas list');
            atlasList = ' ';
        end

        if contains(leadFieldFile, 'Talairach')
            iVal = 2;
        end
    else
        atlasList = ' ';
    end
    userdata{2} = atlasList;
    set(fig, 'userdata', userdata);
    set(findobj(fig, 'tag', 'atlaslist'), 'string', atlasList, 'value', iVal);
    return;
end             
        
if ~isfield(EEG(1).dipfit, 'sourcemodel') || isempty(EEG(1).dipfit.sourcemodel)
    strDipfit   = 'Use DIPFIT leadfield matrix (none present right now)';
    defaultFile = '';
else
    strDipfit   = 'Use pre-calculated DIPFIT Leadfield matrix';
    defaultFile = EEG(1).dipfit.sourcemodel.file;
end
strComputeShort = { 'LCMV' 'LCMVFieldtrip' 'eLoreta' 'eLoretaFieldtrip' };

if nargin < 2
    
    options = {};
    if EEG(end).trials == 1
        if EEG(end).srate > 128
            res = questdlg2( [ 'This function will resample data at 100 Hz, delete all events,' 10 ...
                               'and extract 2-sec data segments. Do you want to proceed?' ], 'Warning ROI connect', 'Cancel', 'Yes', 'Yes');
            if strcmpi(res, 'Cancel'), return; end
            if strcmpi(res, 'yes')
                options = { options{:} 'resample' 100 };
            end
        else
            res = questdlg2( [ 'This function will delete all events and extract 2-sec' 10 ...
                             'data segments. Do you want to proceed?' ], 'Warning ROI connect', 'Cancel', 'Yes', 'Yes');
            if strcmpi(res, 'Cancel'), return; end
            if strcmpi(res, 'yes')
                options = { options{:} };
            end
        end
    elseif EEG(end).srate > 128
        res = questdlg2( 'This function will resample data at 100 Hz. Do you want to proceed?', 'Warning ROI connect', 'Cancel', 'Yes', 'Yes');
        if strcmpi(res, 'Cancel'), return; end
        if strcmpi(res, 'yes')
            options = { options{:} 'resample' 100 };
        end
    end
                     
    leadfield = [];
    leadfield(end+1).label = 'Leadfield matrix: compute using head and source model above (Fieldtrip)';
    leadfield(end+1).label = 'Leadfield matrix: compute using head and source model above (mkfilt_eloreta_v2 - uses less RAM)';
    leadfield(end+1).label = 'Leadfield matrix: Load precomputed leadfield matrix from Fieldtrip, Brainstorm or NFT';
    
    cb_select = 'pop_roi_activity(''cb_select'', gcbf);';
    cb_load   = 'pop_roi_activity(''cb_load'', gcbf);';

    rowg = [0.1 0.5 1 0.2];
    rowg2 = [0.1 1 0.2 0.5];
    strLeadfield = { strDipfit 'Use custom source model aligned to MNI (Brainstorm, Fieldtrip etc...)' };
    strCompute   = { 'Compute distributed source solution using ROIconnect LCMV' ...
                     'Compute distributed source solution using Fieldtrip LCMV' ...
                     'Compute distributed source solution using ROIconnect eLoreta' ...
                     'Compute distributed source solution using Fieldtrip eLoreta'  };
    uigeom = { 1 1 rowg rowg 1 1 rowg 1 1 rowg2 };
    uilist = { { 'style' 'text' 'string' 'Head and source model parameters' 'fontweight' 'bold'} ...
        { 'style' 'popupmenu' 'string' strLeadfield 'tag' 'leadfieldselect' 'callback' cb_select }  ...
        {} { 'style' 'text' 'string' 'Source model file:'  } { 'style' 'edit' 'string' defaultFile 'tag' 'leadfield'   'enable'  'off'   } { 'style' 'pushbutton' 'string' '...' 'tag' 'but' 'callback' cb_load }  ...
        {} { 'style' 'text' 'string' 'Choose ROI atlas:' } { 'style' 'popupmenu' 'string' 'xxxx' 'tag' 'atlaslist' } {} ...
        {} ...
        { 'style' 'popupmenu' 'string' strCompute 'tag' 'model' }  ...
        {} { 'style' 'text' 'string' 'Model parameters:' } { 'style' 'edit' 'string' '0.05'  'tag' 'modelparams' } {}  ...
        {} ...
        { 'style' 'text' 'string' 'Other parameters' 'fontweight' 'bold'} ...
        {} { 'style' 'text' 'string' 'Number of dimensions per ROI:'  } { 'style' 'edit' 'string' '3' 'tag' 'pca' } {}  ...
        };
    
    [result,usrdat,~,out] = inputgui('geometry', uigeom, 'uilist', uilist, 'helpcom', 'pophelp(''pop_roi_activity'')', ...
        'title', 'Compute ROI activity', 'userdata', { defaultFile '' }, 'eval', 'pop_roi_activity(''cb_select'', gcf);');
    if isempty(result), return, end

    % 
    if out.leadfieldselect == 2
         options = { options{:} 'leadfield' EEG(1).dipfit.leadfield }; % file provided by user
         % otherwise use EEG.dipfit.sourcemodel
    end
    try
        modelParams = eval( [ '{' out.modelparams '}' ] );
    catch
        error('Model parameters badly formated');
    end
    options = { options{:} ...
        'model' strComputeShort{ out.model } ...
        'modelparams'  modelParams ...
        'atlas' usrdat{2}{out.atlaslist} ...
        'nPCA'  str2num(out.pca) ...
        };
else
    options = varargin;
end

% process multiple datasets
% -------------------------
if length(EEG) > 1
    if nargin < 2
        [ EEG, com ] = eeg_eval( 'pop_roi_activity', EEG, 'warning', 'on', 'params', options );
    else
        [ EEG, com ] = eeg_eval( 'pop_roi_activity', EEG, 'params', options );
    end
    return;
end

%    'export2icamatrix' 'string' {'on', 'off'}   'off';
[g, moreargs] = finputcheck(options, { ...
    'leadfield'       { 'string' 'struct' } { { }  {} }      '';
    'model'           'string'              strComputeShort  'LCMV';
    'modelparams'     'cell'                {}               {};
    'atlas'           'string'              {}               '';
    'resample'        'real'                {}               100;
    'regepochs'       'string'              { 'on' 'off'}    'off'; % ignored
    'effectchanges'   'string'              { 'on' 'off'}    'off';
    'nPCA'            'real'                {}               3;
    'epochlen'        'real'                {}               2;
    'epochrecur'      'real'                {}               2;
    'fooof'           'string'              { 'on' 'off'}    'off';
    'freqresolution'   'integer'            {}               0;
    'lowmemory'         'string'             { 'on' 'off'}    'off';              
    'fooof_frange'     ''                   {}               [1 30]}, 'pop_roi_activity', 'ignore');
if ischar(g), error(g); end

EEGOUT = EEG;
if ~isempty(g.resample) && ~isequal(EEG.srate, g.resample)
    EEGOUT = pop_resample(EEGOUT, g.resample);
end
if EEGOUT.trials == 1
    EEGOUT = eeg_regepochs(EEGOUT, g.epochrecur, [0 g.epochlen]);
end

chansel = {};
if isempty(g.leadfield)
    g.leadfield = EEG.dipfit.sourcemodel;
    chansel     = EEG.dipfit.sourcemodel.label;
end
if isstruct(g.leadfield) && isfield(g.leadfield, 'file')
    sourceModelFile = g.leadfield.file;
    sourceModel2MNI = g.leadfield.coordtransform;
elseif ~isempty(g.leadfield)
    sourceModelFile = g.leadfield;
    sourceModel2MNI = [];
else
    sourceModelFile = g.leadfield.file;
    sourceModel2MNI = g.leadfield.coordtransform;
    g.leadfield = EEG.dipfit.sourcemodel;
end    

EEGOUT = roi_activity(EEGOUT, 'leadfield', g.leadfield, 'headmodel', EEG.dipfit.hdmfile, ...
    'model', g.model, 'modelparams', g.modelparams, 'sourcemodel', sourceModelFile, ...
    'sourcemodel2mni', sourceModel2MNI, 'nPCA', g.nPCA,'fooof', g.fooof, 'fooof_frange', g.fooof_frange, ...
    'freqresolution', g.freqresolution, 'sourcemodelatlas', g.atlas, 'lowmemory', g.lowmemory, 'chansel', chansel, moreargs{:});

if strcmpi(g.effectchanges, 'on')
    EEG = EEGOUT;
else
    EEG.roi = EEGOUT.roi;
end

if nargout > 1
    for iOption = 1:2:length(options)
        if strcmpi(options{iOption}, 'leadfield') && isequal(options{iOption+1}, EEG.dipfit.sourcemodel)
            options{iOption+1} = '''EEG.dipfit.sourcemodel''';
            break;
        end
    end
    com = sprintf( 'EEG = pop_roi_activity(EEG, %s);', vararg2str( options ));
    com = regexprep(com, '''''''', '');
end

