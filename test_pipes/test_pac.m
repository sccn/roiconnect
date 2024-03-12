% Test PAC implementation for different inputs: single frequencies and frequency bands.
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

%% Test bispectrum for single frequency inputs
low = 10;
high = 50;

fcomb.low = low;
fcomb.high = high;


EEG1 = pop_roi_connect(EEG, 'methods', {'PAC', 'MIM', 'COH'}, 'fcomb', fcomb); % test all 3 connectivity functions (data2spwctrgc, data2strgcmim, roi_pac)
tic
EEG2 = pop_roi_connect(EEG, 'methods', {'PAC'}, 'fcomb', fcomb, 'bs_outopts', 1, 'conn_stats', 'on', 'nshuf', 3, 'poolsize', 12); % compute only b_anti, b_anti_norm
toc
EEG3 = pop_roi_connect(EEG, 'methods', {'PAC'}, 'fcomb', fcomb, 'bs_outopts', 1); % compute only b_anti, b_anti_norm

if ~isequal(squeeze(EEG2.roi.PAC.b_orig(:, :, 1)), EEG3.roi.PAC.b_orig) 
    error 'The first shuffle in the surrogate connectivity array is not the true matrix.'
end

if ~isequal(squeeze(EEG2.roi.PAC.b_anti(:, :, 1)), EEG3.roi.PAC.b_anti) 
    error 'The first shuffle in the surrogate connectivity array is not the true matrix.'
end


%% Test bispectrum for frequency band inputs
low = [4 8];
high = [48 50];

fcomb.low = low;
fcomb.high = high;

tic
EEG4 = pop_roi_connect(EEG, 'methods', {'PAC', 'MIM', 'COH'}, 'fcomb', fcomb); % test all 3 connectivity functions (data2spwctrgc, data2strgcmim, roi_pac)toc
toc
tic
EEG5 = pop_roi_connect(EEG, 'methods', {'PAC'}, 'fcomb', fcomb, 'conn_stats', 'off', 'nshuf', 3, 'bs_outopts', 1); 
toc

if ~isequal(squeeze(EEG5.roi.PAC.b_anti(:, :, 1)), EEG4.roi.PAC.b_anti) 
    error 'The first shuffle in the surrogate connectivity array is not the true matrix.'
end

%% Test PAC plotting
% Test for single frequency inputs
pop_roi_connectplot(EEG3, 'measure', 'pac', 'plotmatrix', 'on');
pop_roi_connectplot(EEG3, 'measure', 'pac_anti', 'plotmatrix', 'on');

% Provoke errors by plotting bispectral tensors that do not exist
pop_roi_connectplot(EEG3, 'measure', 'pac_anti', 'plotmatrix', 'on'); % bs_outopts 4 means only original bispectra are computed, so cannot plot anti
pop_roi_connectplot(EEG3, 'measure', 'pac', 'plotmatrix', 'on'); % bs_outopts 5 means only antisymm. bispectra are computed, so cannot plot normal bispectrum

% Test for frequency bands
pop_roi_connectplot(EEG5, 'measure', 'pac', 'plotmatrix', 'on');
pop_roi_connectplot(EEG5, 'measure', 'pac_anti', 'plotmatrix', 'on');

% plot MIM and COH as a sanity check
pop_roi_connectplot(EEG1, 'measure', 'mim', 'plotmatrix', 'on');
pop_roi_connectplot(EEG1, 'measure', 'coh', 'plotmatrix', 'on');

