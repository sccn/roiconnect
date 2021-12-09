% example of processing ds001787

% start EEGLAB
clear
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

% Install required plugins
plugin_askinstall('bids-matlab-tools', 'pop_importbids', 1);
plugin_askinstall('ICLabel', 'pop_icflag', 1);
plugin_askinstall('iirfilt', 'iirfilt', 1);
plugin_askinstall('PICARD', 'picard', 1);
plugin_askinstall('BIOSIG', 'pop_biosig', 1);
pop_editoptions('option_parallel', 0);

% call BIDS tool BIDS
filepath        = pwd; % use current folder
[STUDY, ALLEEG] = pop_importbids('ds001787', 'studyName','Meditation','outputdir','derivatives', ...
    'eventtype','value','bidsevent','off','bidschanloc','off');
sess2 = [ALLEEG.session] ~= 1;
ALLEEG(sess2) = [];
STUDY.datasetinfo(sess2) = [];
for iDat = 1:length(ALLEEG)
    STUDY.datasetinfo(iDat).index = iDat;
end
ALLEEG = eeg_checkset(ALLEEG, 'loaddata');
CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = [1:length(EEG)];
pop_savestudy(STUDY, ALLEEG, 'filename', 'meditation.study');

%% Preprocessing
% remove channels with no coordinates
EEG = pop_select( EEG, 'nochannel',{'EXG1','EXG2','EXG3','EXG4','EXG5','EXG6','EXG7','EXG8','GSR1','GSR2','Erg1','Erg2','Resp','Plet','Temp'});
for iEEG = 1:length(EEG)
    EEG(iEEG).chanlocs = struct('labels', { 'Fp1' 'AF7' 'AF3' 'F1' 'F3' 'F5' 'F7' 'FT7' 'FC5' 'FC3' 'FC1' 'C1' 'C3' 'C5' 'T7' 'TP7' 'CP5' 'CP3' 'CP1' 'P1' 'P3' 'P5' 'P7' 'P9' 'PO7' 'PO3' 'O1' 'Iz' 'Oz' 'POz' 'Pz' 'CPz' 'Fpz' 'Fp2' 'AF8' 'AF4' 'AFz' 'Fz' 'F2' 'F4' 'F6' 'F8' 'FT8' 'FC6' 'FC4' 'FC2' 'FCz' 'Cz' 'C2' 'C4' 'C6' 'T8' 'TP8' 'CP6' 'CP4' 'CP2' 'P2' 'P4' 'P6' 'P8' 'P10' 'PO8' 'PO4' 'O2' });
end
EEG = eeg_checkset(EEG);
eeglabp = fileparts(which('eeglab.m'));
EEG = pop_chanedit(EEG, 'lookup',fullfile(eeglabp, 'plugins','dipfit','standard_BEM','elec','standard_1005.elc'));
chanlocs = EEG(1).chanlocs;

% Rereference using average reference
EEG = pop_reref( EEG,[]);

% Remove bad channels and bad data
EEG = pop_iirfilt( EEG, 0, 45, [], 0, 0);
EEG = pop_clean_rawdata( EEG,'FlatlineCriterion',5,'ChannelCriterion',0.85,...
    'LineNoiseCriterion',2.5,'Highpass',[0.75 1.25] ,...
    'BurstCriterion',20,'WindowCriterion',0.25,'BurstRejection','off',...
    'Distance','Euclidian','WindowCriterionTolerances',[-Inf 7]);
EEG = pop_rmdat( EEG, {'condition 128' 'condition 1' 'condition 2' 'condition 4' 'condition 8' 'condition 254' '128' '1' '2' '4' '8' '254' },[-1 3] ,1); % remove data around events

% Run ICA and flag artifactual components using IClabel
EEG = pop_runica(EEG, 'icatype','picard', 'pca', -1, 'maxiter', 500);
EEG = pop_iclabel(EEG,'default');
EEG = pop_icflag(EEG,[NaN NaN;0.8 1;0.8 1;NaN NaN;NaN NaN;NaN NaN;NaN NaN]);
EEG = pop_subcomp(EEG, [], 0); % remove pre-flagged bad components
ALLEEG = EEG;

% resample extract epochs
for iEEG = 1:length(ALLEEG)
    EEG = ALLEEG(iEEG);
    EEG = eeg_checkset(EEG, 'loaddata');
    EEG = pop_resample(EEG, 100);
    EEG = eeg_regepochs(EEG, 2, [0 2]);

    % SHOULD WE INTERPOLATE MISSING CHANNELS????????
    EEG = pop_interp(EEG, chanlocs);
    EEG = pop_reref( EEG,[]);

    % Align with head model
    EEG = pop_dipfit_settings( EEG, 'hdmfile',fullfile(eeglabp, 'plugins','dipfit','standard_BEM','standard_vol.mat'), ...
        'coordformat','MNI','mrifile',fullfile(eeglabp, 'plugins','dipfit','standard_BEM','standard_mri.mat'),...
        'chanfile',fullfile(eeglabp, 'plugins','dipfit','standard_BEM','elec', 'standard_1005.elc'),...
        'coord_transform',[0 0 0 0 0 -1.5708 1 1 1] ,'chansel',[1:EEG.nbchan] );

    % compute lead field
    EEG = pop_leadfield(EEG, 'sourcemodel',fullfile(eeglabp,'functions','supportfiles','head_modelColin27_5003_Standard-10-5-Cap339.mat'), ...
        'sourcemodel2mni',[0 -24 -45 0 0 -1.5708 1000 1000 1000] ,'downsample',1);

    % compute activity and connectivity
    EEG = pop_roi_activity(EEG, 'leadfield',EEG.dipfit.sourcemodel,'model','LCMV','modelparams',{0.05},'atlas','LORETA-Talairach-BAs','nPCA',3);
    EEG = pop_roi_connect(EEG, 'methods', { 'CS' });
    
    % resave
    EEG.saved = 'no';
    pop_saveset(EEG, 'savemode', 'resave');
end
