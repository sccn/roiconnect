% Test bispectral time-delay estimation.
%% Run pipeline
clear
eeglab

eeglabp = fileparts(which('eeglab.m'));
EEG = pop_loadset('filename','eeglab_data_epochs_ica.set','filepath',fullfile(eeglabp, 'sample_data/'));
EEG = pop_resample( EEG, 100);
EEG = pop_epoch( EEG, { }, [-0.5 1.5], 'newname', 'EEG Data epochs epochs', 'epochinfo', 'yes');
EEG = pop_select( EEG, 'trial',1:30);
% [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG);
eeglab redraw;

EEG = pop_dipfit_settings( EEG, 'hdmfile',fullfile(eeglabp, 'plugins','dipfit','standard_BEM','standard_vol.mat'), ...
    'coordformat','MNI','mrifile',fullfile(eeglabp, 'plugins','dipfit','standard_BEM','standard_mri.mat'),...
    'chanfile',fullfile(eeglabp, 'plugins','dipfit','standard_BEM','elec', 'standard_1005.elc'),...
    'coord_transform',[0.83215 -15.6287 2.4114 0.081214 0.00093739 -1.5732 1.1742 1.0601 1.1485] );

% EEG = pop_leadfield(EEG, 'sourcemodel',fullfile(eeglabp,'plugins','dipfit','LORETA-Talairach-BAs.mat'), ...
%     'sourcemodel2mni',[0 -24 -45 0 0 -1.5708 1000 1000 1000] ,'downsample',1);

EEG = pop_leadfield(EEG, 'sourcemodel',fullfile(eeglabp,'functions','supportfiles','head_modelColin27_5003_Standard-10-5-Cap339.mat'), ...
    'sourcemodel2mni',[0 -24 -45 0 0 -1.5708 1000 1000 1000] ,'downsample',1);

EEG = pop_roi_activity(EEG, 'leadfield',EEG.dipfit.sourcemodel,'model','LCMV','modelparams',{0.05},'atlas','LORETA-Talairach-BAs','nPCA', 3, 'chansel', EEG.dipfit.chansel);

%% Compute iCOH to extract a promising region combination that shows connectivity
EEG = pop_roi_connect(EEG, 'methods', {'iCOH'});
pop_roi_connectplot(EEG, 'measure', 'iCOH', 'plotcortex', 'off', 'plotmatrix', 'on', 'plotbutterfly', 'on'); % look at broadband iCOH
pop_roi_connectplot(EEG, 'measure', 'iCOH', 'freqrange', [11 13], 'plotcortex', 'off', 'plotmatrix', 'on', 'plotbutterfly', 'on'); 

%% Test TDE on broadband
% find regions from "visual analysis" above
roi1 = 'lingual L'; 
roi2 = 'lateraloccipital L';
labels = replace_underscores(get_labels(EEG));
roi1_idx = find(ismember(labels, roi1));
roi2_idx = find(ismember(labels, roi2));

EEG1 = pop_roi_connect(EEG, 'methods', {'TDE'}, 'tde_regions', [roi1_idx roi2_idx], 'tde_method', 1);

%% Test TDE on frequency bands
EEG2 = pop_roi_connect(EEG, 'methods', {'TDE'}, 'tde_regions', [roi1_idx roi2_idx], 'tde_method', 1, 'tde_freqbands', [11 13]);

%% Plotting
% broadband
pop_roi_connectplot(EEG1, 'measure', 'tde');
pop_roi_connectplot(EEG1, 'measure', 'tde_anti');

% frequency band
pop_roi_connectplot(EEG2, 'measure', 'tde');
pop_roi_connectplot(EEG2, 'measure', 'tde_anti');