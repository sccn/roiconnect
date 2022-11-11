% Requires EEGLAB, ROIconnect, and SIFT installed (SIFT is used to simulate the sources)
% Simulate 2 sources in the temporal cortex
%
% Arnaud Delorme, November 2022

clear
eeglab

%% create empty dataset and create leadfield matrix
eeglabp = fileparts(which('eeglab'));
sourceModelPath = fullfile(eeglabp, 'functions','supportfiles','head_modelColin27_5003_Standard-10-5-Cap339.mat');
tmpSource   = load('-mat', sourceModelPath);
chans =  { 'AF7' 'FP1' 'Afz' 'FP2' 'AF8' 'F7' 'F3' 'F1' 'Fz' 'F2' 'F4' 'F8' 'FT7' 'FC5' 'FC3' 'FC1' 'FCz' 'FC2' 'FC4' 'FC6' 'FT8' 'T7' 'C5' 'C3' 'C1' 'Cz' 'C2' 'C4' 'C6' 'T8' 'TP7' 'CP5' 'CP3' 'CP1' 'CPz' 'CP2' 'CP4' 'CP6' 'TP8' 'P7' 'P5' 'P3' 'P1' 'Pz' 'P2' 'P4' 'P6' 'P8' 'P9' 'PO7' 'PO5' 'PO3' 'Poz' 'PO4' 'PO6' 'PO8' 'P10' 'PO9' 'PO9h' 'O1' 'Oz' 'O2' 'PO10h' 'PO10' };
scalpData   = zeros(length(chans), 200, 30);
EEG = pop_importdata('dataformat','array','nbchan',0,'data',scalpData,'srate',100,'pnts',0,'xmin',0);
EEG.chanlocs = struct('labels', chans);
dipfitdefs;
EEG=pop_chanedit(EEG, 'lookup',template_models(2).chanfile);
EEG = pop_dipfit_settings( EEG, 'hdmfile',template_models(2).hdmfile,'coordformat','MNI',...
    'mrifile',template_models(2).mrifile,'chanfile',template_models(2).chanfile,'coord_transform',[0 0 0 0 0 -1.5708 1 1 1] ,'chansel',[1:EEG.nbchan] );
EEG = pop_leadfield(EEG, 'sourcemodel',fullfile(eeglabp, 'functions/supportfiles/head_modelColin27_5003_Standard-10-5-Cap339.mat'),'sourcemodel2mni',[0 -24 -45 0 0 -1.5708 1000 1000 1000] ,'downsample',1);
leadfield = [ EEG.dipfit.sourcemodel.leadfield{:} ];

%% Simulate VAR model 2 sources coupled using SIFT
%simulations = hlp_getSimExamples;simulations{5}{2} = feval(simulations{5}{2});
simulations{1} = 'Bivariate Coupled Oscillator';
simulations{2} = {'x1(t) = 0.6*x1(t-1) +  0.65*x2(t-2)+ e1(t)' 'x2(t) = 0.5*x2(t-1) + -0.3*x2(t-2) + e2(t)'};
sourceRois = {'temporalpole L' 'temporalpole R' };
[EEGsim] = sim_varmodel('sim', simulations, 'plotData', false, 'simParams', { 'TrialLength', 2, 'NumTrials',30 });

%% project sources onto scalp
if length(sourceRois) ~= size(EEGsim.data,1), error('Wrong number of sources'); end
nvert = size(tmpSource.cortex.vertices,1);
sigma = 10/1000; % same unit as vertices (source spread)
snr   = 20; % signal noise ratio of soucess
P     = zeros(nvert,1);           	% source spatial distribution
J     = zeros(nvert,EEG.pnts*EEG.trials);       % current density estimate
scalpData = zeros(EEG.nbchan,EEG.pnts*EEG.trials);

for iRoi = 1:length(sourceRois)
    % find region
    indRoi = strmatch(sourceRois{iRoi}, tmpSource.atlas.label, 'exact'); 
    x0  = mean(tmpSource.cortex.vertices(tmpSource.atlas.colorTable == indRoi,:)); % location of source
    d   = sqrt(sum((tmpSource.cortex.vertices - ones(nvert,1)*x0).^2,2)); % distance of all vertices to source

    % sum projected current over sources
    P = normpdf(d,0,sigma)*sigma*sqrt(2*pi);  % unit amplitude
    P = repmat(P', [3 1]); P = P(:); % repeat 3 times because 3 orientations
    J = P*EEGsim.data(iRoi,:);
    scalpData = scalpData + leadfield*J;
end
scalpData = scalpData + (std(scalpData(:))./snr).*randn(EEG.nbchan,EEG.pnts*EEG.trials); % add noise
EEG.data  = scalpData;
EEG = eeg_checkset(EEG);

%% Use ROIconnect to find sources
EEG = pop_roi_activity(EEG, 'leadfield',EEG.dipfit.sourcemodel,'model','LCMV','modelparams',{0.05},'atlas','Desikan-Kilianny','nPCA',3);
EEG = pop_roi_connect(EEG, 'morder',2,'naccu',[],'methods',{'CS','COH','MIM'});
pop_roi_connectplot(EEG, 'measure','MIM','freqrange',[],'plotcortex','on','plotcortexparams',{},'plotcortexseedregion',0,'plotmatrix','on','plotpsd','off','plot3d','off','plot3dparams',{'thresholdper',0.8},'region','all','hemisphere','all');
pop_roi_connectplot(EEG, 'measure','MIM','freqrange',[],'plotcortex','off','plotcortexparams',{},'plotcortexseedregion',0,'plotmatrix','off','plotpsd','off','plot3d','on','plot3dparams',{'thresholdper',0.98},'region','all','hemisphere','all');
