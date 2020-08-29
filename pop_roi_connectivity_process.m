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
    fig = EEG;
    userdat = get(fig, 'userdata');
    
    if strcmpi(varargin{1}, 'scale')
        res = inputgui('uilist', { {'style', 'text', 'string' ['Volumetric atlases often have million of voxels and ' 10 ...
            'require to be downsampled to a few thousand voxels' 10 ...
            'to be used as source models. Cancel not to downsample.' 10 ...
            'This message box does not apply to surface atlases.' ]} ...
            {'style', 'text', 'string' 'Downsample (with 4, 1M -> ~15k)' } ...
            {'style', 'edit', 'string' '4' }}, ...
            'geomvert', [3 1], 'geometry', {1 [3 1]});
        if ~isempty(res)
            set(findobj(fig, 'tag', 'down'), 'string', sprintf('Scale/%s', res{1}));
        else
            set(findobj(fig, 'tag', 'down'), 'string', '');
        end
        if ~isempty(res)
            userdat{4} = str2num(res{1});
        else
            userdat{4} = 1;
        end
        set(fig, 'userdata', userdat);
        
    elseif strcmpi(varargin{1}, 'select1') % headmodel
        usrdat = userdat{1}(get(findobj(fig, 'tag', 'selection1'), 'value'));
        set(findobj(fig, 'userdata', 'headmodel'), 'enable', usrdat.enable);
        set(findobj(fig, 'tag', 'strfile1')  , 'string', usrdat.file);
        set(findobj(fig, 'tag', 'transform1'), 'string', num2str(usrdat.align));
        
    elseif strcmpi(varargin{1}, 'select2') % leadfield
        usrdat = userdat{2}(get(findobj(gcf, 'tag', 'selection2'), 'value'));
        set(findobj(gcf, 'userdata', 'leadfield'), 'enable', usrdat.enable);
        set(findobj(gcf, 'tag', 'strfile2')  , 'string', usrdat.file);
        set(findobj(gcf, 'tag', 'selection1')  , 'value', usrdat.headmodel);
        set(findobj(gcf, 'tag', 'selection3')  , 'value', usrdat.sourcemodel);
        pop_roi_connectivity_process(fig, 'select1');
        pop_roi_connectivity_process(fig, 'select3');
        set(findobj(gcf, 'userdata', 'headmodel')  , 'enable', usrdat.modelenable);
        set(findobj(gcf, 'userdata', 'sourcemodel'), 'enable', usrdat.modelenable);
        set(findobj(gcf, 'tag', 'selection1'),     'enable', usrdat.modelenable);
        set(findobj(gcf, 'tag', 'selection3'),     'enable', usrdat.modelenable);
        set(findobj(gcf, 'tag', 'down'),           'visible', 'off');
        
    elseif strcmpi(varargin{1}, 'select3') % atlas
        usrdat = userdat{3}(get(findobj(gcf, 'tag', 'selection3'), 'value'));
        set(findobj(gcf, 'tag', 'push3')     , 'enable', usrdat.enable);
        set(findobj(gcf, 'tag', 'strfile3')  , 'string', usrdat.file, 'enable', usrdat.enable);
        set(findobj(gcf, 'tag', 'transform3'), 'string', num2str(usrdat.align), 'enable', 'on'); % usrdat.enable );
        set(findobj(gcf, 'tag', 'atlas')     , 'string', usrdat.atlasliststr, 'value', usrdat.atlasind, 'enable', 'on' );
        if usrdat.scale>1, set(findobj(gcf, 'tag', 'down'), 'string', sprintf('Scale/%d', usrdat.scale), 'enable', 'on' ); end
        userdat{4} = usrdat.scale;
        set(gcf, 'userdata', userdat);
        
    elseif strcmpi(varargin{1}, 'load1') % headmodel
        [tmpfilename, tmpfilepath] = uigetfile('*', 'Select a text file');
        if tmpfilename(1) ~=0, set(findobj('parent', gcbf, 'tag', 'strfile1'), 'string', fullfile(tmpfilepath,tmpfilename)); end
        
    elseif strcmpi(varargin{1}, 'load2') % leadfield
        [tmpfilename, tmpfilepath] = uigetfile('*', 'Select a text file');
        if tmpfilename(1) ~=0
            fullFileName = fullfile(tmpfilepath,tmpfilename);
            set(findobj('parent', gcbf, 'tag', 'strfile2'), 'string', fullFileName); 
            [~,~,ext] = fileparts(fullFileName);

            % Brain strom files have a link to the source file
            if strcmpi(ext, '.mat')
                warning off
                tmpData = load('-mat', fullFileName, 'SurfaceFile', 'LFM');
                warning on
                if isfield(tmpData, 'SurfaceFile')
                    disp('Brainstorm leadfield matrix detected');
                    
                    % Brain strom files have a link to the source file
                    p = fileparts(fileparts(fileparts(fullFileName)));
                    try
                        surfaceFile = fullfile(p, 'anat', tmpData.SurfaceFile);
                        f = load('-mat', surfaceFile);
                        set(findobj(fig, 'tag', 'strfile3'), 'string', surfaceFile);
                        [ atlasliststr, atlaslist] = getatlaslist(surfaceFile)
                        set(findobj(fig, 'tag', 'atlas')   , 'string', atlasliststr, 'value', 1);
                    catch
                        disp('Could not find Brainstorm anatomical file');
                    end
                    
                elseif isfield(tmpData, 'LFM')
                    disp('NFT leadfield matrix detected');
                    
                    res = questdlg2( [ 'A NFT leadfield matrix file was detected. You should now' 10 'select the associated dipole file to be able to plot the solution' ],  'NFT Dipole File', 'Cancel', 'Continue', 'Continue');
                    % NFT files require an addition .dip source file
                    if strcmpi(res, 'continue')
                        [tmpfilename, tmpfilepath] = uigetfile('*.dip', 'Select the source model file');
                        if tmpfilename(1) ~=0
                            fullFileName = fullfile(tmpfilepath,tmpfilename);
                            set(findobj(fig, 'tag', 'strfile3'), 'string', fullFileName);
                        end
                    end
                end
            end
            
        end
        
    elseif strcmpi(varargin{1}, 'load3') % atlas
        [tmpfilename, tmpfilepath] = uigetfile('*', 'Select a text file');
        if tmpfilename(1) ~=0, set(findobj('parent', gcbf, 'tag', 'strfile3'), 'string', fullfile(tmpfilepath,tmpfilename)); end
        pop_roi_connectivity_process(gcbf, 'scale');
        
    elseif strcmpi(varargin{1}, 'selectcoreg1')
        EEG = userdat{5};
        tmpmodel = get( findobj(gcbf, 'tag', 'strfile1'), 'string');
        tmptransf = get( findobj(gcbf, 'tag', 'transform1'), 'string');
        coregister(EEG.chanlocs, [], 'mesh', tmpmodel,'transform', str2num(tmptransf), 'manual', 'show', 'showlabels1', 'on', 'title', 'Use DIPFIT settings to adjust co-registration');
        
    elseif strcmpi(varargin{1}, 'selectcoreg2')
        plot3dmeshalign(get( findobj(fig, 'tag', 'strfile1'), 'string'), get( findobj(fig, 'tag', 'strfile3'), 'string'), str2num(get( findobj(fig, 'tag', 'transform3'), 'string')));
    end
    return
    
end

if nargin < 2
    
    dipfitOK = false;
    if isfield(EEG.dipfit, 'coordformat')
        dipfitOK = strcmpi(EEG.dipfit.coordformat, 'MNI');
    end
    
    headmodel = [];
    if ~dipfitOK
        res = questdlg2( strvcat('You may use the DIPFIT MNI head model for ROI', ...
            'connectivity analysis. However, you need to go back', ...
            'to the DIPFIT settings to align it with your montage.', ...
            'To continue, you must have a custom Leadfield matrix.'), 'Use DIPFIT Leadfield matrix', 'Continue', 'Go back', 'Go back');
        if strcmpi(res, 'Go back'), return; end
    else
        headmodel(1).label = 'Headmodel: Use DIPFIT current model and conductivities';
        headmodel(1).file  = EEG.dipfit.hdmfile;
        headmodel(1).align = EEG.dipfit.coord_transform;
        headmodel(1).enable = 'off';
    end
    headmodel(end+1).label = 'Headmodel: DIPFIT compatible head model & conductivities in MNI space';
    headmodel(end).file  = '';
    headmodel(end).align = [];
    headmodel(end).enable = 'on';
    
    leadfield = [];
    leadfield(end+1).label = 'Leadfield matrix: compute using head model and source model above (Fieldtrip)';
    leadfield(end).enable  = 'off';
    leadfield(end).file  = '';
    leadfield(end).headmodel = 1;
    leadfield(end).sourcemodel = 1;
    leadfield(end).modelenable = 'on';

    leadfield(end+1).label = 'Leadfield matrix: compute using head model and source model above (mkfilt_eloreta_v2)';
    leadfield(end).enable  = 'off';
    leadfield(end).file  = '';
    leadfield(end).headmodel = 1;
    leadfield(end).sourcemodel = 1;
    leadfield(end).modelenable = 'on';

    leadfield(end+1).label = 'Leadfield matrix: Load precomputed leadfield matrix from Fieldtrip, Brainstorm or NFT';
    leadfield(end).enable  = 'on';
    leadfield(end).file  = '';
    leadfield(end).headmodel = 2;
    leadfield(end).sourcemodel = 4;
    leadfield(end).modelenable = 'off';
    
    p  = fileparts(which('eeglab.m'));
    roi(1).label = 'Source model ROI: Use Desikan-Kilianny in Colin27 template';
    roi(1).file  = fullfile( p, 'functions', 'supportfiles', 'head_modelColin27_5003_Standard-10-5-Cap339.mat');
    roi(1).align = [0 -24 -45 0 0 -1.5707963 1000 1000 1000];
    roi(1).enable = 'off';
    roi(1).scale  = 1;
    roi(1).atlasliststr = { 'Desikan-Kiliany (68 ROIs)' };
    roi(1).atlaslist    = { 'Desikan-Kiliany' };
    roi(1).atlasind  = 1;
    
    p  = fileparts(which('pop_roi_connectivity_process.m'));
    roi(2).label = 'Source model ROI: Use Desikan-Kilianny in ICBM152 template (Brainstrom)';
    roi(2).file  = fullfile(p, 'tess_cortex_mid_low_2000V.mat');
    roi(2).align = [0 -24 -45 0 0 -1.5707963000 1000 1000 1000];
    roi(2).enable = 'off';
    roi(2).scale  = 1;
    [ roi(2).atlasliststr, roi(2).atlaslist] = getatlaslist(roi(2).file);
    roi(2).atlasind  = 2;
    
    roi(3).label = 'Source model ROI: LORETA-KEY';
    roi(3).file  = fullfile(p, 'LORETA-Talairach-BAs.mat');
    roi(3).align = [];
    roi(3).enable = 'off';
    roi(3).scale  = 1;
    roi(3).atlasliststr = { 'LORETA-Talairach-BAs (44 x 2 ROIs)' };
    roi(3).atlaslist    = { 'LORETA-Talairach-BAs' };
    roi(3).atlasind  = 1;
        
    p  = fileparts(which('ft_defaults.m'));
    roi(4).label = 'Source model ROI: AFNI TTatlas+tlrc (Fieldtrip)';
    roi(4).file  = fullfile(p, 'template','atlas','afni','TTatlas+tlrc.HEAD');
    roi(4).align = [ ];
    roi(4).enable = 'off';
    roi(4).scale  = 4;
    roi(4).atlasliststr = { '' };
    roi(4).atlaslist    = { '' };
    roi(4).atlasind  = 1;

    roi(5).label = 'Source model ROI: Custom atlas';
    roi(5).file  = '';
    roi(5).align = [];
    roi(5).enable = 'on';
    roi(5).scale  = 4;
    roi(5).atlasliststr = { '' };
    roi(5).atlaslist    = { '' };
    roi(5).atlasind     = 1;
    
    cb_select1 = 'pop_roi_connectivity_process(gcf, ''select1'');';
    cb_select2 = 'pop_roi_connectivity_process(gcf, ''select2'');';
    cb_select3 = 'pop_roi_connectivity_process(gcf, ''select3'');';
    cb_load1   = 'pop_roi_connectivity_process(gcf, ''load1'');';
    cb_load2   = 'pop_roi_connectivity_process(gcf, ''load2'');';
    cb_load3   = 'pop_roi_connectivity_process(gcf, ''load3'');';
    cb_selectcoreg1 = 'pop_roi_connectivity_process(gcf, ''selectcoreg1'');';
    cb_selectcoreg2 = 'pop_roi_connectivity_process(gcf, ''selectcoreg2'');';
    
    rowg = [0.1 0.6 1 0.2];
    uigeom = { 1 1 rowg rowg 1 rowg rowg [0.1 0.6 0.9 0.3] 1 rowg 1 [0.5 1 0.35 0.5] [0.5 1 0.35 0.5] [1] [0.2 1 1.5] };
    uilist = { { 'style' 'text' 'string' 'Region Of Interest (ROI) connectivity analysis' 'fontweight' 'bold'} ...
        { 'style' 'popupmenu' 'string' { headmodel.label } 'tag' 'selection1' 'callback' cb_select1 }  ...
        {} { 'style' 'text' 'string' 'File name' 'userdata', 'headmodel'             } { 'style' 'edit' 'string' 'xxxxxxxxxxxxxxxxxxxx' 'tag' 'strfile1'   'enable'  'off' 'userdata', 'headmodel'   } { 'style' 'pushbutton' 'string' '...' 'userdata', 'headmodel' 'tag' 'push1' 'callback' cb_load1 }  ...
        {} { 'style' 'text' 'string' 'Elec to MNI align' 'userdata', 'headmodel' } { 'style' 'edit' 'string' 'xxxxxxxxxxxxxxxxxxxx' 'tag' 'transform1' 'enable'  'off' 'userdata', 'headmodel' } { 'style' 'pushbutton' 'string' '...' 'callback' cb_selectcoreg1 } ...
        ...
        { 'style' 'popupmenu' 'string' { roi.label }  'tag' 'selection3' 'callback' cb_select3 } ...
        {} { 'style' 'text' 'string' 'File name'}                { 'style' 'edit' 'string' 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx' 'tag' 'strfile3'   'userdata', 'sourcemodel' } { 'style' 'pushbutton' 'string' '...'  'userdata', 'sourcemodel' 'tag' 'push3' 'callback' cb_load3 }  ...
        {} { 'style' 'text' 'string' 'File to headmodel align' } { 'style' 'edit' 'string' 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx' 'tag' 'transform3' 'userdata', 'sourcemodel' } { 'style' 'pushbutton' 'string' '...'  'userdata', 'sourcemodel' 'callback' cb_selectcoreg2 }  ...
        {} { 'style' 'text' 'string' 'Use this Atlas/ROI'      } { 'style' 'popupmenu' 'string' 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx' 'tag' 'atlas'                           } { 'style' 'text' 'string' 'Scale/2'    'userdata', 'sourcemodel' 'tag' 'down' }  ...
        { 'style' 'popupmenu' 'string' { leadfield.label } 'tag' 'selection2' 'callback' cb_select2 }  ...
        {} { 'style' 'text' 'string' 'File name'  'userdata', 'leadfield' } { 'style' 'edit' 'string' 'xxxxxxxxxxxxxxxxxxxx' 'tag' 'strfile2' 'enable'  'off' 'userdata', 'leadfield' } { 'style' 'pushbutton' 'string' '...' 'tag' 'push2' 'callback' cb_load2 'userdata', 'leadfield' }  ...
        ...
        {} ...
        {} { 'style' 'text' 'string' 'Autoregressive model order' } { 'style' 'edit' 'string' '20' 'tag' 'morder' } { } ...
        {} { 'style' 'text' 'string' 'Bootstrap if any (n)' } { 'style' 'edit' 'string' '' 'tag' 'naccu' } { } ...
        {} ...
        {} { 'style' 'checkbox' 'string' 'Compute TRGC' 'tag' 'trgc' 'value' 1 } ...
        { 'style' 'checkbox' 'string' 'Compute cross-spectrum' 'tag' 'crossspec' 'value' 1 } ...
        };
    [result,usrdat,~,out] = inputgui('geometry', uigeom, 'uilist', uilist, 'helpcom', 'pophelp(''pop_loadbv'')', ...
        'title', 'Load a Brain Vision Data Exchange format dataset', 'userdata', {headmodel leadfield roi [] EEG}, 'eval', [cb_select1 cb_select2 cb_select3 'set(findobj(gcf, ''tag'', ''down''), ''string'', '''');' ]);
    if isempty(result), return, end
    
    if isempty(usrdat{3}), usrdat{3} = 1; end
    options = {};
    options = { 'headmodel' out.strfile1 ...
        'elec2mni' str2num(out.transform1) ...
        'leadfield' out.strfile2 ...
        'sourcemodel' out.strfile3 ...
        'sourcemodel2mni' str2num(out.transform3) ...
        'sourcemodelatlas' roi(out.selection3).atlaslist{out.atlas} ...
        'sourceanalysis' fastif(out.selection2 == 1, 'fieldtrip', 'roiconnect') ...
        'downsample' usrdat{4} ...
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
    'leadfield'       'string'  { }             '';
    'elec2mni'        'real'    { }             [];
    'downsample'      'integer'  { }             4; % volume only
    'sourcemodel'     'string'  { }             '';
    'sourcemodel2mni' 'real'    { }             [] }, 'pop_roi_connectivity_process', 'ignore');
if ischar(g), error(g); end

if isempty(g.leadfield)
    % Source model
    headmodel = load('-mat', g.headmodel);
    EEG.dipfit.coord_transform = g.elec2mni;
    dataPre = eeglab2fieldtrip(EEG, 'preprocessing', 'dipfit'); % does the transformation
    ftPath = fileparts(which('ft_defaults'));

    % Prepare the liedfield matrix
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
    elseif strcmpi(ext, '.mat') % && isfield(g.sourcemodel, 'tri')
        sourcemodelOri = transform_move_inward(g.sourcemodel, g.headmodel,g.sourcemodel2mni);
    end

    cfg      = [];
    cfg.elec = dataPre.elec;
    %     cfg.grid    = sourcemodelOri;   % source points
    if isfield(headmodel, 'vol')
        cfg.headmodel = headmodel.vol;   % volume conduction model
    else
        cfg.headmodel = headmodel;   % volume conduction model
    end
    cfg.sourcemodel.inside = ones(size(sourcemodelOri.pos,1),1) > 0;
    cfg.sourcemodel.pos    = sourcemodelOri.pos;
    if isfield(sourcemodelOri, 'tri')
        cfg.sourcemodel.tri    = sourcemodelOri.tri;
    end
    cfg.singleshell.batchsize = 5000; % speeds up the computation
    leadfield = ft_prepare_leadfield(cfg);
else
    leadfield = g.leadfield;
    if isstr(leadfield)
        leadfield = load('-mat', leadfield);
    end
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

EEG = roi_connectivity_process(EEG, 'leadfield', leadfield, 'headmodel', g.headmodel, 'sourcemodel', g.sourcemodel, 'downsample', g.downsample, 'sourcemodel2mni', g.sourcemodel2mni, moreargs{:});

if nargout > 1
    com = sprintf( 'EEG = pop_roi_connectivity_process(EEG, %s);', vararg2str( options ));
end