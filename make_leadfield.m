clear

EEG = pop_loadset(['analysis_output/sub_NDARAA948VFH/preprocessing/data/prep_interp.set']);
EEG = pop_dipfit_settings( EEG, 'hdmfile','/data/matlab/eeglab/plugins/dipfit3.3/standard_BEM/standard_vol.mat','coordformat','MNI','mrifile','/data/matlab/eeglab/plugins/dipfit/standard_BEM/standard_mri.mat','chanfile','/data/matlab/eeglab/plugins/dipfit/standard_BEM/elec/standard_1005.elc','coord_transform',[0 -18 0 0 0 -1.571 11.3431 11.3431 11.3431] ,'chansel',[1:128] );
dataPre = eeglab2fieldtrip(EEG, 'preprocessing', 'dipfit');

headmodel = load('-mat', EEG.dipfit.hdmfile);
ftPath = fileparts(which('ft_defaults'));

sourcemodelOriOld = ft_read_headshape(fullfile(ftPath, 'template', 'sourcemodel', 'cortex_20484.surf.gii'));
p  = fileparts(which('eeglab.m'));
hm = load('-mat', fullfile( p, 'functions', 'supportfiles', 'head_modelColin27_5003_Standard-10-5-Cap339.mat'));
clear sourcemodelOri2;
tf = traditionaldipfit([0.0000000000 -26.6046230000 -46.0000000000 0.1234625600 0.0000000000 -1.5707963000 1000.0000000000 1000.0000000000 1000.0000000000]);
sourcemodelOri.pos      = tf*[hm.cortex.vertices ones(size(hm.cortex.vertices,1),1)]';
sourcemodelOri.pos      = sourcemodelOri.pos';
sourcemodelOri.pos(:,4) = [];
sourcemodelOri.tri = hm.cortex.faces;
sourcemodelOri.unit = 'mm';
% sourcemodelOri = sourcemodelOriOld;
%     figure; plot3dmesh(sourcemodelOriOld.pos, sourcemodelOriOld.tri);
%     plot3dmesh(sourcemodelOri.pos, sourcemodelOri.tri);
%     title('Superposed meshes');

cfg         = [];
cfg.elec            = dataPre.elec;
%     cfg.grid    = sourcemodelOri;   % source points
cfg.headmodel = headmodel.vol;   % volume conduction model
cfg.sourcemodel.inside = ones(size(sourcemodelOri.pos,1),1) > 0;
cfg.sourcemodel.pos    = sourcemodelOri.pos;
cfg.sourcemodel.tri    = sourcemodelOri.tri;
cfg.singleshell.batchsize = 5000; % speeds up the computation
sourcemodel = ft_prepare_leadfield(cfg);
save('-mat', 'leadfield_tmp.mat', '-struct', 'sourcemodel');

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
p  = fileparts(which('eeglab.m'));
cortexfile = fullfile( p, 'functions', 'supportfiles', 'head_modelColin27_5003_Standard-10-5-Cap339.mat');
roi_connectivity_process(EEG, 'leadfield', 'leadfield_tmp.mat', 'cortexfile', cortexfile);
