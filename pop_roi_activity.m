% pop_roi_activity - call roi_activity to compute activities of ROIs using
%                    eLoreta
%
% Usage:
%  EEG = pop_roi_activity(EEG, 'key', 'val', ...);
%
% Inputs:
%  EEG - EEGLAB dataset
%
% Required inputs:
%  'headmodel'   - [string] head model file in MNI space
%  'sourcemodel' - [string] source model file
%
% Optional inputs:
%  'elec2mni'    - [9x float] homogeneous transformation matrix to convert
%                  electrode locations to MNI space.
%  'sourcemodel2mni' - [9x float] homogeneous transformation matrix to convert
%                  sourcemodel to MNI space.
%
% Output:
%  EEG - EEGLAB dataset with field 'roi' containing connectivity info.
%
% Note: Optional inputs to roi_activity() are also accepted.
%
% Author: Arnaud Delorme, UCSD, 2019
%
% Example
%   p = fileparts(which('eeglab')); % path
%   EEG = pop_roi_activity(EEG, 'headmodel', ...
%   EEG.dipfit.hdmfile, 'elec2mni', EEG.dipfit.coord_transform, ...
%   'sourcemodel', fullfile(p, 'functions', 'supportfiles', ...
%   'head_modelColin27_5003_Standard-10-5-Cap339.mat'), 'sourcemodel2mni', ...
%   [0 -26.6046230000 -46 0.1234625600 0 -1.5707963000 1000 1000 1000]);
%
% Use pop_roi_connect(EEG) to compute conectivity

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
    else
        atlasList = ' ';
    end
    userdata{2} = atlasList;
    set(fig, 'userdata', userdata);
    set(findobj(fig, 'tag', 'atlaslist'), 'string', atlasList);
    return;
end             
        
if ~isfield(EEG.dipfit, 'sourcemodel') || isempty(EEG.dipfit.sourcemodel)
    strDipfit   = 'Use DIPFIT leadfield matrix (none present right now)';
    defaultFile = '';
else
    strDipfit   = 'Use pre-calculated DIPFIT Leadfield matrix';
    defaultFile = EEG.dipfit.sourcemodel.file;
end
strComputeShort = { 'LCMV' 'LCMVFieldtrip' 'eLoreta' 'eLoretaFieldtrip' };

if nargin < 2
    
    options = {};
    if EEG.trials == 1
        if EEG.srate > 128
            res = questdlg2( [ 'This function is optimized to process 2-sec data epochs' 10 ...
                'at a sampling rate of about 100 Hz. Do you want to resample the data and' 10 ...
                'extract 2-sec data segments? (make sure your dataset is saved)' ], 'Warning ROI connect', 'Cancel', 'No', 'Yes', 'Yes');
            if strcmpi(res, 'cancel'), return; end
            if strcmpi(res, 'yes')
                options = { options{:} 'resample' 'on' 'regepochs' 'on' };
            end
        else
            res = questdlg2( [ 'This function is optimized to process 2-sec data epochs.' 10 ...
                'Do you want to resample the data extract 2-sec data segments?' 10 ...
                '(make sure your dataset is saved)' ], 'Warning ROI connect', 'Cancel', 'No', 'Yes', 'Yes');
            if strcmpi(res, 'cancel'), return; end
            if strcmpi(res, 'yes')
                options = { options{:} 'regepochs' 'on' };
            end
        end
    elseif EEG.srate > 128
        res = questdlg2( [ 'This function is optimized to process data epochs' 10 ...
            'at a sampling rate of about 100 Hz. Do you want to resample the data?' 10 ...
            '(make sure your dataset is saved)' ], 'Warning ROI connect', 'Cancel', 'No', 'Yes', 'Yes');
        if strcmpi(res, 'cancel'), return; end
        if strcmpi(res, 'yes')
            options = { options{:} 'resample' 'on' };
        end
    end
                     
    leadfield = [];
    leadfield(end+1).label = 'Leadfield matrix: compute using head and source model above (Fieldtrip)';
    leadfield(end+1).label = 'Leadfield matrix: compute using head and source model above (mkfilt_eloreta_v2 - uses less RAM)';
    leadfield(end+1).label = 'Leadfield matrix: Load precomputed leadfield matrix from Fieldtrip, Brainstorm or NFT';
    
    cb_select = 'pop_roi_activity(''cb_select'', gcbf);';
    cb_load   = 'pop_roi_activity(''cb_load'', gcbf);';

    rowg = [0.1 0.5 1 0.2];
    strLeadfield = { strDipfit 'Use custom source model aligned to MNI (Brainstorm, Fieldtrip etc...)' };
    strCompute   = { 'Compute distributed source solution using ROIconnect LCMV' ...
                     'Compute distributed source solution using Fieldtrip LCMV' ...
                     'Compute distributed source solution using ROIconnect eLoreta' ...
                     'Compute distributed source solution using Fieldtrip eLoreta'  };
    uigeom = { 1 1 rowg rowg 1 1 rowg [0.1 0.5 0.2 1] };
    uilist = { { 'style' 'text' 'string' 'Region Of Interest (ROI) connectivity analysis' 'fontweight' 'bold'} ...
        { 'style' 'popupmenu' 'string' strLeadfield 'tag' 'leadfieldselect' 'callback' cb_select }  ...
        {} { 'style' 'text' 'string' 'Source model file:'  } { 'style' 'edit' 'string' defaultFile 'tag' 'leadfield'   'enable'  'off'   } { 'style' 'pushbutton' 'string' '...' 'tag' 'but' 'callback' cb_load }  ...
        {} { 'style' 'text' 'string' 'Choose ROI atlas:' } { 'style' 'popupmenu' 'string' 'xxxx' 'tag' 'atlaslist' } {} ...
        {} ...
        { 'style' 'popupmenu' 'string' strCompute 'tag' 'model' }  ...
        {} { 'style' 'text' 'string' 'Model parameters:' } { 'style' 'edit' 'string' '0.05'  'tag' 'modelparams' } {}  ...
        {} { 'style' 'text' 'string' 'Dimention per ROI:'  } { 'style' 'edit' 'string' '3' 'tag' 'pca' } {}  ...
        };
    
    [result,usrdat,~,out] = inputgui('geometry', uigeom, 'uilist', uilist, 'helpcom', 'pophelp(''pop_roi_activity'')', ...
        'title', 'Compute ROI activity', 'userdata', { defaultFile '' }, 'eval', 'pop_roi_activity(''cb_select'', gcf);');
    if isempty(result), return, end

    % 
    if out.leadfieldselect == 1
         options = { options{:} 'leadfield' EEG.dipfit.sourcemodel };
    else
         options = { options{:} 'leadfield' EEG.dipfit.leadfield };
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
    % check that the dipfit settings are the same
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
    'model'           'string'              strComputeShort 'LCMV';
    'modelparams'     'cell'                {}               {};
    'atlas'           'string'              {}               '';
    'resample'        'string'              { 'on' 'off'}    'off';
    'regepochs'       'string'              { 'on' 'off'}    'off';
    'nPCA'            'real'                {}               3 }, 'pop_roi_activity', 'ignore');
if ischar(g), error(g); end

if strcmpi(g.resample, 'on')
    EEG = pop_resample(EEG, 100);
end
if strcmpi(g.regepochs, 'on')
    EEG = eeg_regepochs(EEG, 2, [0 2]);
end

if isstruct(g.leadfield) && isfield(g.leadfield, 'file')
    sourceModelFile = g.leadfield.file;
    sourceModel2MNI = g.leadfield.coordtransform;
else
    sourceModelFile = g.leadfield;
    sourceModel2MNI = [];
end    
%modelParams;

EEG = roi_activity(EEG, 'leadfield', g.leadfield, 'headmodel', EEG.dipfit.hdmfile, ...
    'model', g.model, 'modelparams', g.modelparams, 'sourcemodel', sourceModelFile, ...
    'sourcemodel2mni', sourceModel2MNI, 'nPCA', g.nPCA, ...
    'sourcemodelatlas', g.atlas, moreargs{:});

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

