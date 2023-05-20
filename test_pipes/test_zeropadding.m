% Test zero-padding after choosing a (custom) frequency resolution (in number 
% of frequencies) in roi_activitiy and roi_connect.
%% Run pipeline
clear
eeglab

eeglabp = fileparts(which('eeglab.m'));
EEG = pop_loadset('filename','eeglab_data_epochs_ica.set','filepath',fullfile(eeglabp, 'sample_data/'));
srate = 100; % can be tested with different srates but they are always resampled to 100 Hz
EEG = pop_resample( EEG, srate);
EEG = pop_epoch( EEG, { }, [-0.5 1.5], 'newname', 'EEG Data epochs epochs', 'epochinfo', 'yes');
EEG = pop_select( EEG, 'trial',1:30);
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG);
eeglab redraw;

EEG = pop_dipfit_settings( EEG, 'hdmfile',fullfile(eeglabp, 'plugins','dipfit','standard_BEM','standard_vol.mat'), ...
    'coordformat','MNI','mrifile',fullfile(eeglabp, 'plugins','dipfit','standard_BEM','standard_mri.mat'),...
    'chanfile',fullfile(eeglabp, 'plugins','dipfit','standard_BEM','elec', 'standard_1005.elc'),...
    'coord_transform',[0.83215 -15.6287 2.4114 0.081214 0.00093739 -1.5732 1.1742 1.0601 1.1485] ,'chansel',[1:32] );

EEG = pop_leadfield(EEG, 'sourcemodel',fullfile(eeglabp,'plugins','dipfit','LORETA-Talairach-BAs.mat'), ...
    'sourcemodel2mni',[0 -24 -45 0 0 -1.5708 1000 1000 1000] ,'downsample',1);

%% Test zero padding in roi_activity and plot results
% sanity check with default frequency resolution (freqresolution = 0)
EEG = pop_roi_activity(EEG, 'leadfield',EEG.dipfit.sourcemodel,'model','LCMV','modelparams',{0.05},'atlas','LORETA-Talairach-BAs','nPCA',3);

% subplot(2,1,1)
% for j = 1:size(EEG.roi.source_roi_power,2)
%     plot(EEG.roi.freqs,EEG.roi.source_roi_power(:,j))
%     hold on;
% end
% original_nfreq = length(EEG.roi.freqs);
% title("Number of frequencies: " + int2str(original_nfreq) + " (Default)")
% xlabel('Frequency [Hz]')

% specify desired frequency resolution (in number of desired frequencies)
desired_nfreq = 400;
EEG = pop_roi_activity(EEG, 'leadfield',EEG.dipfit.sourcemodel,'model','LCMV','modelparams',{0.05},'atlas','LORETA-Talairach-BAs','nPCA',3,'freqresolution',desired_nfreq);

% subplot(2,1,2)
% for j = 1:size(EEG.roi.source_roi_power,2)
%     plot(EEG.roi.freqs,EEG.roi.source_roi_power(:,j))
%     hold on;
% end
% title("Number of frequencies: " + int2str(desired_nfreq + 1))
% xlabel('Frequency [Hz]')

%% Perform zero padding in roi_connect and plot results
% test data2sctrgcmim
EEG = pop_roi_connect(EEG, 'methods', {'MIM'}); % default (without padding)
EEG = pop_roi_connect(EEG, 'methods', {'MIM'}, 'freqresolution', desired_nfreq); % with padding

% test data2spwctrgc
EEG = pop_roi_connect(EEG, 'methods', { 'CS', 'COH'}, 'snippet', 'on', 'snip_length', 60, 'fcsave_format', 'mean_snips'); % also test data2spwctrgc
EEG = pop_roi_connect(EEG, 'methods', { 'CS', 'COH'}, 'snippet', 'on', 'snip_length', 60, 'fcsave_format', 'mean_snips', 'freqresolution', desired_nfreq);

