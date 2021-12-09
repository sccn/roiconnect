%% directories

dropbox_dir = '/Volumes/data/Users/haufe/Dropbox/Zurich_Networks/';

data_dir = [dropbox_dir 'EEG_data/data_for_test_retest/'];

brainstorm_db_dir = [data_dir 'brainstorm_db/Zurich_test_retest/'];

bs_result_dir = [data_dir 'processed_bs/'];

fig_folder = '../../figures/source_imaging/';
sbj = 'A23';
mkdir([fig_folder sbj])


%% settings

target_fs = 100;

% frequency resolution using 100 bins (actually 101, from 0 to 50 in steps
% of 0.5 Hz)
fres = target_fs;

% length of epochs in sec
epochLeng = 2;

% 2 sec epochs for spectral/connectivity analysis
winlength = epochLeng*target_fs;

% model order for the spectral/connectivity analysis
morder = 20;

eloreta_reg = 0.05;
lcmv_reg = 0.05;

% either an integer for the exact number of components
% or a percentage of variance explained
nPCA = 3;

% colormap
load cm17;

% plot smooth cortical surfaces
smooth_cortex = 0.35;

% some parameters for connectivity estimation
lpFilter =   45;       % low-pass filter cut-off
bsFilter =   [47 53];       % additional notch filter 
filterOrder = 2;    % Butterworth filter order
dsRatio =  5;       % downsampling rate
epochLeng = 2;   % epoch length in seconds

% from the MVGC toolbox, compute frequencies in Hz
frqs = sfreqs(fres, target_fs);

% frequency bin indices for some EEG rhythms
bands = {[4 8], [8 13], [13 28]};
bandnames = {'theta', 'alpha', 'beta'};
for iba = 1:length(bands)
  frq_inds{iba} = find(frqs >= bands{iba}(1) & frqs < bands{iba}(2));
end

% number of bootstrap samples
nbootstrap = [];

%% load data
load([bs_result_dir 'bs_results']);

cortex.Vertices = cortex.Vertices*1000;

load([data_dir 'A23/ggiip_A23_RES_EEG.mat']);
load([brainstorm_db_dir 'data/A23/@default_study/channel.mat']);

%% zero-phase filtering
[b_low, a_low] = butter(5, lpFilter/(EEG.srate/2), 'low');
[b_notch, a_notch] = butter(filterOrder, bsFilter/(EEG.srate/2),'stop');
a_all = poly([roots(a_low);roots(a_notch)]);
b_all = conv(b_low,b_notch);
% fvtool(b_all, a_all)
EEG.data = filtfilt(b_all, a_all, double(EEG.data)')';
    
%% plain downsampling to 100 Hz. I prefer this to resampling to 100 Hz
%% to avoid all the interpolation steps involved in up and downsampling
EEG.data = EEG.data(:, 1:dsRatio:end);
EEG.srate = EEG.srate/dsRatio;
EEG.pnts    = size(EEG.data,2);
EEG.xmax    = EEG.xmin + (EEG.pnts-1)/EEG.srate;
EEG.times   = linspace(EEG.xmin*1000, EEG.xmax*1000, EEG.pnts);
EEG = eeg_checkset( EEG );


[Nchan, Ntp]=size(EEG.data);
EEG = pop_select( EEG,'point',[1 (Ntp - mod(Ntp,EEG.srate*epochLeng))] );
Nepoch = EEG.pnts / (EEG.srate*epochLeng);
EEG = eeg_checkset(EEG);
%EEG.data = EEG.data(:, 1:(Ntp - mod(Ntp,fsNew*8)) );
for ievent = 1: Nepoch
    EEG.event(ievent).type = num2str(epochLeng);
    EEG.event(ievent).latency = (ievent-1)*(EEG.srate*epochLeng)+1;
    EEG.event(ievent).duration = epochLeng;
end
EEG = eeg_checkset( EEG );
EEG = pop_epoch( EEG, {  num2str(epochLeng)  }, [0  epochLeng], 'epochinfo', 'yes');

% EEG.data = double(EEG.data);

%% source imaging

% check
assert(isequal({Channel.Name}, {EEG.chanlocs.labels}))

% rereference data and leadfield to common average
H = eye(EEG.nbchan) - ones(EEG.nbchan) ./ EEG.nbchan;

EEG.data = reshape(H*EEG.data(:, :), EEG.nbchan, EEG.pnts, EEG.trials);
leadfield = reshape(H*leadfield(:, :), EEG.nbchan, nvox, 3);

% number of ROIs in the Desikan-Kiliany Atlas
nROI = length(cortex.Atlas(3).Scouts);

% ROI labels
labels = {cortex.Atlas(3).Scouts.Label};

ind_cortex = [];
ind_roi = {};
nvoxel_per_roi = [];
for iROI = 1:nROI
  ind_roi{iROI} = cortex.Atlas(3).Scouts(iROI).Vertices;
  ind_cortex = cat(1, ind_cortex, ind_roi{iROI});
  [~, ind_roi_cortex{iROI}, ~] = intersect(ind_cortex, ind_roi{iROI});
  nvoxel_per_roi(iROI) = length(ind_roi_cortex{iROI});
end
nvox = length(ind_cortex);

leadfield = leadfield(:, ind_cortex, :);

% eLORETA projection kernel
% P = mkfilt_eloreta_v2(leadfield, eloreta_reg);

% LCMV projection kernel, works better than eloreta in our simulations
C = cov(EEG.data(:, :)');
alpha = lcmv_reg*trace(C)/length(C);
Cr = C + alpha*eye(EEG.nbchan);
[~, P] = lcmv(Cr, leadfield, struct('alpha', 0, 'onedim', 0));

% project to source space
% the units are nA*m
source_voxel_data = reshape(EEG.data(:, :)'*P(:, :), EEG.pnts*EEG.trials, nvox, 3);
source_voxel_data = 10^3*source_voxel_data;

% new: sensor_space cross-spectrum
% used to compute full ROI wise spectra later
conn = data2spwctrgc(reshape(EEG.data, [], EEG.pnts, EEG.trials), fres, morder, 0, [], [], {'CS'});
CS_sensor = conn.CS;

% keep only the first nPCA strongest components for each ROI
% where nPCA is either fixed or determined by a percentage of explained variance (e.g. 90%)
source_roi_data = [];
source_roi_power = [];
source_roi_power_norm = [];
source_roi_power_total = [];
source_roi_power_total_norm = [];
varex = ones(nROI, EEG.nbchan);
nPCAs = [];
for iROI = 1:nROI
  data_ = source_voxel_data(:, ind_roi_cortex{iROI}, :);
  va_ = var(data_(:, :));
  source_roi_power(iROI) = sum(va_)';
  % power normalized per area
  source_roi_power_norm(iROI) = source_roi_power(iROI)/length(ind_roi_cortex{iROI});
  
  % optional z-scoring, this makes the PCA independent of the power in each
  % voxel, and favors to find components that are consistently expressed in
  % many voxels rather than only in a few voxels with strong power (which
  % may leak from a neighboring region) 
  % removed because it did not yield an advantage in our simulations 
  %   data_(:, :) = zscore(data_(:, :))*sqrt(mean(va_));

  % new way to compute ROI spectral power: base it on all voxels/directions
  % instead of inly the first nPCA PCA components
  P_ = P(:, ind_roi_cortex{iROI}, :);
  P_ = P_(:, :);
  CS_ROI = [];
  for ifreq = 1:size(CS_sensor, 1)
    CS_ROI(ifreq, :, :) = P_'*squeeze(CS_sensor(ifreq, :, :))*P_;
  end
  source_roi_power_total(:, iROI) = sum(cs2psd(CS_ROI), 2);
  source_roi_power_total_norm(:, iROI) = source_roi_power_total(:, iROI)/length(ind_roi_cortex{iROI}); 
 
  % PCA/SVD
  [data_, S_] = svd(data_(:, :), 'econ');
  
  % variance explained
  vx_ = cumsum(diag(S_).^2)./sum(diag(S_).^2);
  invx = 1:min(length(vx_), EEG.nbchan);
  varex(iROI, invx) = vx_(invx);
  
  if nPCA < 1
    nPCAs(iROI) = min(find(varex(iROI, invx) > nPCA)); 
  else
    nPCAs(iROI) = nPCA;
  end
  
  % keep nPCA components with correct unit and scaling
  source_roi_data = cat(2, source_roi_data, data_(:, 1:nPCAs(iROI))*S_(1:nPCAs(iROI), 1:nPCAs(iROI)));
end

% final time series data with nPCA components per ROI
% the dimensions are nROI, nPCA, pnts, trials
% units are still in nA*m

% version with nPCA components
source_roi_data = permute(reshape(source_roi_data, EEG.pnts, EEG.trials, []), [3 1 2]);

% slightly modified version:
% scale up each ROI as if the nPCA component explain 100% of the variance
% of that ROI
% units still in nA*m
% source_roi_data_corr = source_roi_data .* repmat(sqrt(1./varex(:, nPCA)), 1, nPCA, EEG.pnts, EEG.trials);


beg_inds = cumsum([1 nPCAs(1:end-1)]);
end_inds = cumsum([nPCAs]);

for iroi = 1:nROI
  PCA_inds{iroi} = beg_inds(iroi):end_inds(iroi);
end

%% ROI analysis using nPCA components per ROI

% test TRGC between ROIs (that is, pairs of nPCA-dimensional spaces)
% compute indices
% new: 
inds = {}; ninds = 0;
for iroi = 1:nROI
  for jroi = (iroi+1):nROI
    inds{ninds+1} = {PCA_inds{iroi}, PCA_inds{jroi}};    
%     inds{ninds+2} = {PCA_inds{jroi}, PCA_inds{iroi}};  
    ninds = ninds + 1;
  end
end


conn_uni = data2spwctrgc(source_roi_data, fres, morder, 0, nbootstrap, [], {'CS'});


% compute TRGC
% note, this is not yet the net score i->j - j->i , 
% only the difference forward-backward


% TRGC_ = data2sctrgc(source_roi_data, fres, morder, 0, nbootstrap, [], inds);
% new FC measures: MIC and MIM (Ewald and Nolte) are multivariate undirected (symmetric) FC
% metrics. They represent the first eigenvalue and the sum of eigenvalues
% of a maximization of iCOH. 
conn_mult = data2sctrgcmim(source_roi_data, fres, morder, 0, nbootstrap, [], inds, {'TRGC', 'GC', 'MIM', 'MIC'});


%% TRGC plotting

% calculation of net TRGC scores (i->j minus j->i), recommended procedure
% TRGCnet = TRGC_(:, 1:2:end)-TRGC_(:, 2:2:end);
% new way to compute net scores
TRGCnet = conn_mult.TRGC(:, :, 1) - conn_mult.TRGC(:, :, 2);

% create a ROI x ROI connectivity matrix, if needed
% TRGCmat(f, ii, jj) is net TRGC from jj to ii
TRGC = [];
iinds = 0;
for iroi = 1:(nROI)
  for jroi = (iroi+1):(nROI)
    iinds = iinds + 1;
    TRGC(:, iroi, jroi) = -TRGCnet(:, iinds);
    TRGC(:, jroi, iroi) = TRGCnet(:, iinds);
  end
end

% plot overall TRGC spectrum, sanity check for numerical instabilities
figure; plot(frqs, TRGC(:, :))

% TRGC per frequency band
TRGCband = [];
for iba = 1:length(bands)
  TRGCband(iba, :, :) = sum(TRGC(frq_inds{iba}, :, :));
end

% plot alpha-band net TRGC for each region
roi_ = sum(TRGCband(2, :, :), 3);
allplots_cortex_BS(cortex_highres, roi_, [-max(abs(roi_)) max(abs(roi_))], cm17, 'TRGC', smooth_cortex, ...
  [fig_folder sbj '/TRGC_net_alpha']);

figure; imagesc(squeeze(TRGCband(2, :, :)))

% use left inferior parietal region as sender and plot TRGC to all other
% areas
% positive means seed sends more than it receives
roi_ = squeeze(TRGCband(2, :, 15));
allplots_cortex_BS(cortex_highres, roi_, [-max(abs(roi_)) max(abs(roi_))], cm17, 'TRGC', smooth_cortex, ...
  [fig_folder sbj '/TRGC_net_alpha_left_inferior_parietal_sender']);


%% iCOH plotting

% note: as an alternative to iCOH, the MIM/MIC scores could be brought into
% matrix form and be plotted.

% compute (absolute value of) the imaginary part of coherence
iCOH_ = abs(imag(cs2coh(conn_uni.CS)));

clear iCOH
for iroi = 1:nROI
  for jroi = 1:nROI
    iCOH(:, iroi, jroi) = mean(mean(iCOH_(:, PCA_inds{iroi}, PCA_inds{jroi}), 2), 3);
  end
end


% iCOH per frequency band
iCOHband = [];
for iba = 1:length(bands)
  iCOHband(iba, :, :) = sum(iCOH(frq_inds{iba}, :, :));
end

%plot alpha-band net iCOH for each region
roi_ = mean(iCOHband(2, :, :), 3);
allplots_cortex_BS(cortex_highres, roi_, [min(abs(roi_)) max(abs(roi_))], cm17a, '|iCOH|', smooth_cortex, ...
  [fig_folder sbj '/iCOH_net_alpha']);

figure; imagesc(squeeze(iCOHband(2, :, :)))

% use right paracentral region as seed and plot iCOH with all other
% areas
% positive means more interaction
roi_ = squeeze(iCOHband(2, :, 34));
allplots_cortex_BS(cortex_highres, roi_, [min(abs(roi_)) max(abs(roi_))], cm17a, '|iCOH|', smooth_cortex, ...
  [fig_folder sbj '/iCOH_net_alpha_right_paracentral_seed']);


%%

% create a ROI x ROI MIC connectivity matrix
MIC = [];
iinds = 0;
for iroi = 1:(nROI)
  for jroi = (iroi+1):(nROI)
    iinds = iinds + 1;
    MIC(:, iroi, jroi) = conn_mult.MIC(:, iinds);
    MIC(:, jroi, iroi) = conn_mult.MIC(:, iinds);
  end
end

% MIC per frequency band
MICband = [];
for iba = 1:length(bands)
  MICband(iba, :, :) = sum(MIC(frq_inds{iba}, :, :));
end

%plot alpha-band net MIC for each region
roi_ = mean(MICband(2, :, :), 3);
allplots_cortex_BS(cortex_highres, roi_, [min(abs(roi_)) max(abs(roi_))], cm17a, 'MIC', smooth_cortex, ...
  [fig_folder sbj '/MIC_net_alpha']);

% use right paracentral region as seed and plot MIC with all other
% areas
% positive means more interaction
roi_ = squeeze(MICband(2, :, 34));
allplots_cortex_BS(cortex_highres, roi_, [min(abs(roi_)) max(abs(roi_))], cm17a, 'MIC', smooth_cortex, ...
  [fig_folder sbj '/MIC_net_alpha_right_paracentral_seed']);



%% power plotting

% plot regionwise broadband power on smoothed surface
voxel_ = nan*ones(size(cortex.Vertices, 1), 1);
voxel_(ind_cortex) = 10*log10(reshape(sum(sum(source_voxel_data.^2), 3), [], 1));
allplots_cortex_BS(cortex_highres, voxel_(in_normal_to_high), [min(voxel_) max(voxel_)], ...
  cm17a, 'power [dB]', smooth_cortex, [fig_folder sbj '/voxel_broadband_voxel_power_db_smooth']);


% % plot regionwise broadband power
% source_roi_power_norm_dB = 10*log10(source_roi_power_norm);
% allplots_cortex_BS(cortex_highres, source_roi_power_norm_dB, [min(source_roi_power_norm_dB) max(source_roi_power_norm_dB)], ...
%   cm17a, 'power [dB]', 0, ...
%   [fig_folder sbj '/roi_broadband_power_db']);

% plot regionwise broadband power on smoothed surface
source_roi_power_norm_dB = 10*log10(source_roi_power_norm);
allplots_cortex_BS(cortex_highres, source_roi_power_norm_dB, [min(source_roi_power_norm_dB) max(source_roi_power_norm_dB)], ...
  cm17a, 'power [dB]', smooth_cortex, [fig_folder sbj '/roi_broadband_voxel_power_db_smooth']);


clear PS
PS_ = cs2psd(conn_uni.CS);
for iroi = 1:nROI
  PS(:, iroi) = sum(PS_(:, PCA_inds{iroi}), 2);
end

% another way to plot region-wise broadband power 
% todo: need to check for the proper scales after the Fourier transform
roi_ = 10*log10(sum(PS)./nvoxel_per_roi)';
allplots_cortex_BS(cortex_highres, roi_, [min(roi_) max(roi_)], cm17a, 'power [dB]', smooth_cortex, ...
  [fig_folder sbj '/roi_broadband_power_db2']);
      

% power per frequency band
PSband = [];
for iba = 1:length(bands)
  PSband(iba, :) = sum(PS(frq_inds{iba}, :));
end

% 
roi_ = 10*log10(PSdBband(2, :)./nvoxel_per_roi);
allplots_cortex_BS(cortex_highres, roi_, [min(roi_) max(roi_)], cm17a, 'power [dB]', smooth_cortex, ...
  [fig_folder sbj '/power_dB_alpha']);
