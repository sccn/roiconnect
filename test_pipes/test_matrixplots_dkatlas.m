% Test region-to-region FC matrices (MIM) with different parameters. Here, the Desikan-Killiany atlas with 68 ROIs is used as the source model.
%% Run pipeline
clear
eeglab

eeglabp = fileparts(which('eeglab.m'));
EEG = pop_loadset('filename','eeglab_data_epochs_ica.set','filepath',fullfile(eeglabp, 'sample_data/'));
EEG = pop_resample( EEG, 100);
EEG = pop_epoch( EEG, { }, [-0.5 1.5], 'newname', 'EEG Data epochs epochs', 'epochinfo', 'yes');
EEG = pop_select( EEG, 'trial',1:30);
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG);
eeglab redraw;

EEG = pop_dipfit_settings( EEG, 'hdmfile',fullfile(eeglabp, 'plugins','dipfit','standard_BEM','standard_vol.mat'), ...
    'coordformat','MNI','mrifile',fullfile(eeglabp, 'plugins','dipfit','standard_BEM','standard_mri.mat'),...
    'chanfile',fullfile(eeglabp, 'plugins','dipfit','standard_BEM','elec', 'standard_1005.elc'),...
    'coord_transform',[0.83215 -15.6287 2.4114 0.081214 0.00093739 -1.5732 1.1742 1.0601 1.1485] ,'chansel',[1:32] );

EEG = pop_leadfield(EEG, 'sourcemodel',fullfile(eeglabp,'functions','supportfiles','head_modelColin27_5003_Standard-10-5-Cap339.mat'), ...
    'sourcemodel2mni',[0 -24 -45 0 0 -1.5708 1000 1000 1000] ,'downsample',1);

EEG = pop_roi_activity(EEG, 'leadfield',EEG.dipfit.sourcemodel,'model','LCMV','modelparams',{0.05},'atlas','LORETA-Talairach-BAs','nPCA',3);
EEG = pop_roi_connect(EEG, 'methods', {'MIM'});

%% Plot matrix with different parameters
% matrix without any filters
pop_roi_connectplot(EEG, 'measure', 'mim', 'plotcortex', 'off', 'plotmatrix', 'on');

% group by hemispheres
pop_roi_connectplot(EEG, 'measure', 'mim', 'plotcortex', 'off', 'plotmatrix', 'on', 'grouphemispheres', 'on');

% frequency band specified
pop_roi_connectplot(EEG, 'measure', 'mim', 'plotcortex', 'off', 'plotmatrix', 'on', 'freqrange', [4 8]');

% hemisphere specified
pop_roi_connectplot(EEG, 'measure', 'mim', 'plotcortex', 'off', 'plotmatrix', 'on', 'hemisphere', 'right'); 

% cortical region specified
pop_roi_connectplot(EEG, 'measure', 'mim', 'plotcortex', 'off', 'plotmatrix', 'on', 'region', 'occipital');






