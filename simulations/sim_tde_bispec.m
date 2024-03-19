% Simulate artificial two-channel dataset to estimate time delays using 
% antisymmetrized bispectra. The times series consists of a signal and noise 
% component that can be specified.

%% Setup
clear
clc

srate = 100; % Hz
t = 130 * srate; % length of signal in datapoints (130 seconds)
seglen = 2 * srate; % segment length in datapoints (2 seconds)
delay = 50; % delay in time bins

param = 0.7 ;        % parameter used to define distribution. lambda for 
theta = [0.7 0.7];

signal_type = 'exp'; % ["exp", "G", "pink", "pinksq", "spiky"]
noise_type = 'G';  % ["exp", "G", "pink", "pinksq", "spiky"]

snr = 0.8;
beta = 1;

nboot = []; % leave to allow plotting
method = 1; % 1:4, chose which bispectral TDE should be displayed, check documentation in bispectral_TD_est

%% Create data
eeglab
Xraw = mk_series(t+abs(delay), signal_type, 0.7,0);
noiseX = mk_series(t, noise_type, 0.7, 0);
noiseY = mk_series(t, noise_type, 0.7, 0);

% Combine and segment
[Xlong, Ylong] = combine_sn(Xraw, noiseX, noiseY, delay, srate, snr, theta, beta);
X = segment_data(Xlong, seglen);
Y = segment_data(Ylong, seglen);

% Unmixed noise comparison
[XlongU, YlongU] = combine_sn(Xraw, noiseX, noiseY, delay, srate, snr,[0 0], beta);

%% Estimate bispectral TDE
ndat = seglen * 2;
segshift = floor(ndat/2);
epleng = ndat;
fres = srate;
frqs = sfreqs(2 * fres, srate);

[B2_xxx] = squeeze(data2bs_univar(Xlong', 2 * seglen, segshift, epleng, length(frqs)-1));
para_xyx.chancomb = [1, 2, 1]; 
[B2_xyx] = data2bs_univar([Xlong', Ylong'], 2 * seglen, segshift, epleng, length(frqs)-1, para_xyx);
[B2_yyy] = squeeze(data2bs_univar(Ylong', 2 * seglen, segshift, epleng, length(frqs)-1));

% required for antisymmetrization
para_yxx.chancomb = [2, 1, 1]; 
[B2_yxx] = data2bs_univar([Xlong', Ylong'], 2 * seglen, segshift, epleng, length(frqs)-1, para_yxx);

[T, I] = bispectral_TD_est(B2_xxx, B2_xyx, B2_yyy, method, [], 1);
[aT, aI] = bispectral_TD_est(B2_xxx, B2_xyx - B2_yxx, B2_yyy, method, [], 1);

% TDE on frequency bands
band = [20 30];
fmask = true(1, length(frqs));
fmask(frqs < band(1) | frqs > band(2)) = 0;
[T_, I_] = bispectral_TD_est(B2_xxx, B2_xyx, B2_yyy, method, fmask(1:end-1), 1);


%% Plotting
delay_scale = (-seglen+1:seglen-1) / srate;
[peak_val, peak_idx] = max(aT); % extract estimated delay/peak
est_delay = delay_scale(peak_idx);

figure; plot(delay_scale, aT, 'black')
xline(est_delay, '--r')
xlabel('Time (s)')
ylabel('a.u.')
title(sprintf('TDE | Method %d', method))
subtitle("\tau = " + num2str(est_delay) + " (s)")
grid on

%% Check if the estimated delay is in fact the true simulated delay
if ~est_delay == delay
    error('The estimated delay does not match the true simulated delay')
else
    disp('The estimated delay matches the true simulated delay.')
end
