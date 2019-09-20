% pop_roi_connectivity_process - call roi_connectivity_process to compute
%                                connectivity between ROIs
% Usage:
%  EEG = pop_roi_connectivity_process(EEG, 'key', 'val', ...);
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
%  EEG - EEGLAB dataset with field 'roiconnect' containing connectivity info.
%
% Note: Optional inputs to roi_connectivity_process() are also accepted.
%
% Author: Arnaud Delorme, UCSD, 2019
%
% Example
%   p = fileparts(which('eeglab')); % path
%   EEG = pop_roi_connectivity_process(EEG, 'headmodel', ...
%   EEG.dipfit.hdmfile, 'elec2mni', EEG.dipfit.coord_transform, ...
%   'sourcemodel', fullfile(p, 'functions', 'supportfiles', ...
%   'head_modelColin27_5003_Standard-10-5-Cap339.mat'), 'sourcemodel2mni', ...
%   [0 -26.6046230000 -46 0.1234625600 0 -1.5707963000 1000 1000 1000]);
%
% Use pop_roi_connectivity_plot(EEG) to plot the results.

% TO DO - Arno
% - Centralize reading head mesh and Atlas (there might be a function in
% Fieldtrip to do that) ft_read_volume ft_read_mesh
% - Make compatible with all Fieldtrip and FSL Atlases
% - Downsampling of Atlas - check bug submitted to Fieldtrip
% - Plot inside(blue) vs outside(red) voxels for source volume

function [EEG,com] = pop_roi_connectivity_process(EEG, varargin)

com = '';
if nargin < 1
    help pop_roi_connectivity_process;
    return
end

% special callback for custom source models
if ~isstruct(EEG)
    res = inputgui('uilist', { {'style', 'text', 'string' ['Volumetric atlases often have million of voxels and ' 10 ...
                                                          'require to be downsampled to a few thousand voxels' 10 ...
                                                          'to be used as source models. Cancel not to downsample.' 10 ...
                                                          'This message box does not apply to surface atlases.' ]} ...
                                                {'style', 'text', 'string' 'Downsample (with 4, 1M -> ~15k)' } ...
                                                {'style', 'edit', 'string' '4' }}, ...
                                                'geomvert', [3 1], 'geometry', {1 [3 1]});
    if ~isempty(res)
        set(findobj(EEG, 'tag', 'down'), 'string', sprintf('Scale/%s', res{1}));
    else
        set(findobj(EEG, 'tag', 'down'), 'string', '');
    end
    userdata = get(EEG, 'userdata');
    userdata{3} = str2num(res{1});
    set(EEG, 'userdata', userdata);
    return
end

if nargin < 2
    
    dipfitOK = false;
    if isfield(EEG.dipfit, 'coordformat')
        dipfitOK = strcmpi(EEG.dipfit.coordformat, 'MNI');
    end
    
    leadfield = [];
    if ~dipfitOK
        res = questdlg2( strvcat('You may use the DIPFIT MNI head model for ROI', ...
                           'connectivity analysis. However, you need to go back', ...
                           'to the DIPFIT settings to align it with your montage.', ...
                           'To continue, you must have a custom Leadfield matrix.'), 'Use DIPFIT Leadfield matrix', 'Continue', 'Go back', 'Go back');
        if strcmpi(res, 'Go back'), return; end
    else
        leadfield(1).label = 'Leadfield matrix: Use DIPFIT head model in MNI space';
        leadfield(1).file  = EEG.dipfit.hdmfile;
        leadfield(1).align = EEG.dipfit.coord_transform;
        leadfield(1).enable = 'off';
    end
    leadfield(end+1).label = 'Leadfield matrix: Use custom head model in MNI or Brainstrom space';
    leadfield(end).file  = '';
    leadfield(end).align = [];
    leadfield(end).enable = 'on';
    
    p  = fileparts(which('eeglab.m'));
    roi(1).label = 'ROI: Use Desikan-Kilianny in Colin27 template';
    roi(1).file  = fullfile( p, 'functions', 'supportfiles', 'head_modelColin27_5003_Standard-10-5-Cap339.mat');
    roi(1).align = [0 -24 -45 0 0 -1.5707963 1000 1000 1000];
    roi(1).enable = 'off';
    roi(1).atlasliststr = { 'Desikan-Kiliany (68 ROIs)' };
    roi(1).atlaslist    = { 'Desikan-Kiliany' };
    roi(1).atlasind  = 1;
    
    p  = fileparts(which('pop_roi_connectivity_process.m'));
    roi(2).label = 'ROI: Use Desikan-Kilianny in ICBM152 template (Brainstrom)';
    roi(2).file  = fullfile(p, 'tess_cortex_mid_low_2000V.mat');
    roi(2).align = [0 -24 -45 0 0 -1.5707963000 1000 1000 1000];
    roi(2).enable = 'off';
    [ roi(2).atlasliststr, roi(2).atlaslist] = getatlaslist(roi(2).file);
    roi(2).atlasind  = 2;

    roi(3).label = 'ROI: Custom atlas';
    roi(3).file  = '';
    roi(3).align = [];
    roi(3).enable = 'on';
    roi(3).atlasliststr = { '' };
    roi(3).atlaslist    = { '' };
    roi(3).atlasind     = 1;
    
    cb_select1 = [ 'usrdat = get(gcf, ''userdata'');' ...
                  'usrdat = usrdat{1}(get(findobj(gcf, ''tag'', ''selection1''), ''value''));' ...
                  'set(findobj(gcf, ''tag'', ''push1''), ''enable'', usrdat.enable);' ...
                  'set(findobj(gcf, ''tag'', ''strfile1'')  , ''string'', usrdat.file, ''enable'', usrdat.enable);' ...
                  'set(findobj(gcf, ''tag'', ''transform1''), ''string'', num2str(usrdat.align), ''enable'', usrdat.enable );' ...
                  'clear usrdat;' ];
    cb_select2 = [ 'usrdat = get(gcf, ''userdata'');' ...
                  'usrdat = usrdat{2}(get(findobj(gcf, ''tag'', ''selection2''), ''value''));' ...
                  'set(findobj(gcf, ''tag'', ''push2''), ''enable'', usrdat.enable);' ...
                  'set(findobj(gcf, ''tag'', ''strfile2'')  , ''string'', usrdat.file, ''enable'', usrdat.enable);' ...
                  'set(findobj(gcf, ''tag'', ''transform2''), ''string'', num2str(usrdat.align), ''enable'', ''on'');' ... % usrdat.enable );' ...
                  'set(findobj(gcf, ''tag'', ''atlas'')     , ''string'', usrdat.atlasliststr, ''value'', usrdat.atlasind, ''enable'', ''on'' );' ...
                  'clear usrdat;' ];
              
    cb_load1 = [ '[tmpfilename, tmpfilepath] = uigetfile(''*'', ''Select a text file'');' ...
                 'if tmpfilename(1) ~=0, set(findobj(''parent'', gcbf, ''tag'', ''strfile1''), ''string'', fullfile(tmpfilepath,tmpfilename)); end;' ...
                 'clear tmpfilename tmpfilepath;' ];   
             
    cb_load2 = [ '[tmpfilename, tmpfilepath] = uigetfile(''*'', ''Select a text file'');' ...
                 'if tmpfilename(1) ~=0, set(findobj(''parent'', gcbf, ''tag'', ''strfile2''), ''string'', fullfile(tmpfilepath,tmpfilename)); end;' ...
                 'clear tmpfilename tmpfilepath; pop_roi_connectivity_process(gcbf);' ];
             
    cb_selectcoreg1 = [ 'tmpmodel = get( findobj(gcbf, ''tag'', ''strfile1''), ''string'');' ...
                       'tmptransf = get( findobj(gcbf, ''tag'', ''transform1''), ''string'');' ...
                       'coregister(EEG.chanlocs, [], ''mesh'', tmpmodel,''transform'', str2num(tmptransf), ''manual'', ''show'', ''showlabels1'', ''on'', ''title'', ''Use DIPFIT settings to adjust co-registration'');' ...
                       'clear tmpmodel tmptransf;' ];
    cb_selectcoreg2 = [ 'tmpmodel1 = get( findobj(gcbf, ''tag'', ''strfile1''), ''string'');' ...
                        'tmpmodel2 = get( findobj(gcbf, ''tag'', ''strfile2''), ''string'');' ...
                        'tmptransf = str2num(get( findobj(gcbf, ''tag'', ''transform2''), ''string''));' ...
                        'figure; plot3dmeshalign(tmpmodel1);' ...
                        '[tmpmodel2,tmptransf] = transform_move_inward(tmpmodel2, tmpmodel1, tmptransf);' ...
                        'hold on; plot3dmeshalign(tmpmodel2, tmptransf, [1 0 0]);' ...
                        'hlegend = legend({''Head model'' ''ROI source model'' });' ...
                        'set(hlegend, ''position'', [0.7473 0.7935 0.2304 0.0774]);' ...
                        'clear hlegend tmpmodel1 tmpmodel2 tmptransf;' ];
   
    rowg = [0.1 0.6 1 0.2];
    uigeom = { 1 1 rowg rowg 1 rowg rowg [0.1 0.6 0.9 0.3] 1 [0.5 1 0.35 0.5] [0.5 1 0.35 0.5] [1] [0.2 1 1.5] };
    uilist = { { 'style' 'text' 'string' 'Region Of Interest (ROI) connectivity analysis' 'fontweight' 'bold'} ...
              { 'style' 'popupmenu' 'string' { leadfield.label } 'tag' 'selection1' 'callback' cb_select1 }  ...
              {} { 'style' 'text' 'string' 'File name' } { 'style' 'edit' 'string' 'xxxxxxxxxxxxxxxxxxxx' 'tag' 'strfile1' 'enable'  'off' } { 'style' 'pushbutton' 'string' '...' 'tag' 'push1' 'callback' cb_load1 }  ...
              {} { 'style' 'text' 'string' 'Elec to MNI transfrom' } { 'style' 'edit' 'string' 'xxxxxxxxxxxxxxxxxxxx' 'tag' 'transform1' 'enable'  'off'  }  { 'style' 'pushbutton' 'string' '...' 'callback' cb_selectcoreg1 } ...
              { 'style' 'popupmenu' 'string' { roi.label }  'tag' 'selection2' 'callback' cb_select2 } ...
              {} { 'style' 'text' 'string' 'File name' } { 'style' 'edit' 'string' 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx' 'tag' 'strfile2' } { 'style' 'pushbutton' 'string' '...' 'tag' 'push2' 'callback' cb_load2 }  ...
              {} { 'style' 'text' 'string' 'File to MNI transfrom' } { 'style' 'edit' 'string' 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx' 'tag' 'transform2' } { 'style' 'pushbutton' 'string' '...' 'callback' cb_selectcoreg2 }  ...
              {} { 'style' 'text' 'string' 'Use this Atlas/ROI' } { 'style' 'popupmenu' 'string' 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx' 'tag' 'atlas' } { 'style' 'text' 'string' 'Scale/2' 'tag' 'down' }  ...
              {} ...
              {} { 'style' 'text' 'string' 'Model order for AR model' } { 'style' 'edit' 'string' '20' 'tag' 'morder' } { } ...
              {} { 'style' 'text' 'string' 'Bootstrap if any (n)' } { 'style' 'edit' 'string' '' 'tag' 'naccu' } { } ...
              {} ...
              {} { 'style' 'checkbox' 'string' 'Compute TRGC' 'tag' 'trgc' 'value' 1 } ...
              { 'style' 'checkbox' 'string' 'Compute cross-spectrum' 'tag' 'crossspec' 'value' 1 } ...
              };
    [result,usrdat,~,out] = inputgui('geometry', uigeom, 'uilist', uilist, 'helpcom', 'pophelp(''pop_loadbv'')', ...
        'title', 'Load a Brain Vision Data Exchange format dataset', 'userdata', {leadfield roi []}, 'eval', [cb_select1 cb_select2 'set(findobj(gcf, ''tag'', ''down''), ''string'', '''');' ]);
    if isempty(result), return, end

    if ~out.trgc && ~out.crossspec
        error('Nothing to compute')
    end
    if isempty(usrdat{3}), usrdat{3} = 1; end
    options = {};
    options = { 'headmodel' out.strfile1 ...
                'elec2mni' str2num(out.transform1) ...
                'sourcemodel' out.strfile2 ...
                'sourcemodel2mni' str2num(out.transform2) ...
                'sourcemodelatlas' roi(out.selection2).atlaslist{out.atlas} ...
                'downsample' usrdat{3} ...
                'morder' str2num(out.morder) ...
                'naccu' str2num(out.naccu) ...
                'trgc'  fastif(out.trgc, 'on', 'off') ...
                'crossspec' fastif(out.crossspec, 'on', 'off') ...
                }; 
else 
    options = varargin;
end

[g, moreargs] = finputcheck(options, { ...
    'headmodel'       'string'  { }             '';
    'elec2mni'        'real'    { }             [];
    'downsample'      'integer'  { }             4; % volume only
    'sourcemodel'     'string'  { }             '';
    'sourcemodel2mni' 'real'    { }             [] }, 'pop_roi_connectivity_process', 'ignore');
if ischar(g), error(g); end

% Prepare the liedfield matrix
headmodel = load('-mat', g.headmodel);
EEG.dipfit.coord_transform = g.elec2mni;
dataPre = eeglab2fieldtrip(EEG, 'preprocessing', 'dipfit'); % does the transformation
ftPath = fileparts(which('ft_defaults'));

[~,~,ext] = fileparts(g.sourcemodel);
if strcmpi(ext, '.nii')
    atlas = ft_read_atlas(g.sourcemodel);
    mri = sum(atlas.tissue(:,:,:,:),4) > 0;
    [r,c,v] = ind2sub(size(mri),find(mri));
    xyz = [r c v ones(length(r),1)];
    xyz = atlas.transform*xyz';
    if nargin > 1 && ~isempty(transform)
        xyz = traditionaldipfit(transform)*xyz;
    end
    disp('DOWNSAMPLING NOT IMPLEMENTED FOR THIS TYPE OF ATLAS');
elseif strcmpi(ext, '.head')
    [~, sourcemodelOri.pos, ~ ] = load_afni_atlas(g.sourcemodel, g.headmodel, g.sourcemodel2mni, g.downsample);
elseif strcmpi(ext, '.mat')
    sourcemodelOri = transform_move_inward(g.sourcemodel, g.headmodel,g.sourcemodel2mni);
end

if isfield(headmodel, 'vol')
    cfg      = [];
    cfg.elec = dataPre.elec;
    %     cfg.grid    = sourcemodelOri;   % source points
    cfg.headmodel = headmodel.vol;   % volume conduction model
    cfg.sourcemodel.inside = ones(size(sourcemodelOri.pos,1),1) > 0;
    cfg.sourcemodel.pos    = sourcemodelOri.pos;
    if isfield(sourcemodelOri, 'tri')
        cfg.sourcemodel.tri    = sourcemodelOri.tri;
    end
    cfg.singleshell.batchsize = 5000; % speeds up the computation
    sourcemodel = ft_prepare_leadfield(cfg);
else
    sourcemodel = g.headmodel;
end

% remove vertices not modeled (no longer necessary - makes holes in model)
%     indRm = find(sourcemodel.inside == 0);
%     rowRm = [];
%     for ind = 1:length(indRm)
%         sourcemodel.tri(sourcemodel.tri(:,1) == indRm(ind),:) = [];
%         sourcemodel.tri(sourcemodel.tri(:,2) == indRm(ind),:) = [];
%         sourcemodel.tri(sourcemodel.tri(:,3) == indRm(ind),:) = [];
%         sourcemodel.tri(sourcemodel.tri(:) > indRm(ind)) = sourcemodel.tri(sourcemodel.tri(:) > indRm(ind)) - 1;
%     end
%     sourcemodel.pos(indRm,:) = [];
%     sourcemodel.leadfield(indRm) = [];

EEG = roi_connectivity_process(EEG, 'leadfield', sourcemodel, 'headmodel', g.headmodel, 'sourcemodel', g.sourcemodel, 'sourcemodel2mni', g.sourcemodel2mni, moreargs{:});
if nargout > 1
    com = sprintf( 'EEG = pop_roi_connectivity_process(EEG, %s);', vararg2str( options ));
end