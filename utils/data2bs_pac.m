% pac_bispec() - Estimates phase-amplitude coupling (PAC) by computing  
%                bicoherence. Direct relations between PAC and bispectra 
%                are shown in Zandvoort and Nolte, 2021.
%
% Inputs:
%  data   - [array] (nchan x ntimepoints x ntrials) or (nchan x ntimepoints * ntrials) ROI data, requires prior source reconstruction 
%  params - [struct] Struct containing parameters for bispectrum computation
%  
% Outputs:
%   b_orig - ROI x ROI bispectrum
%   b_anti - ROI x ROI antisymmetrized bispectrum
%   b_orig_norm - ROI x ROI bicoherence (normalized by threenorm)
%   b_anti_norm - ROI x ROI antisymmetrized bicoherence (normalized by threenorm)
%
% Copyright (C) Franziska Pellegrini, franziska.pellegrini@charite.de,
%               Tien Dung Nguyen, tien-dung.nguyen@charite.de
%               Zixuan Liu, zixuan.liu@campus.tu-berlin.de

function [b_orig, b_anti, b_orig_norm,b_anti_norm] = data2bs_pac(data, params)

% determine ROIs
nroi = params.nROI;
segleng = params.segleng;
segshift = params.segshift;
epleng = params.epleng;
fcomb = params.fcomb;
fs = params.fs;

fres = fs;
frqs = sfreqs(fres, fs);

% extract all individual frequencies in the selected bands
size_low = size(fcomb.low, 2);
size_high = size(fcomb.high, 2);
mask_inds_low = frqs >= fcomb.low(1) & frqs <= fcomb.low(size_low);
mask_inds_high = frqs >= fcomb.high(1) & frqs <= fcomb.high(size_high);
frqs_low = frqs(mask_inds_low); 
frqs_high = frqs(mask_inds_high);

% determine all frequency combinations
[m, n] = ndgrid(frqs_low, frqs_high);
frqs_combs = [m(:),n(:)]; 
n_combs = size(frqs_combs, 1);
if n_combs > 20
    % according to our test simulations, the computation time scales linearly with the number of frequency pairs times 2, assuming no other ongoing CPU-heavy processes
    time_est = 2 * n_combs; 
    warning('PAC is going to be estimated on %d frequency pair(s). Estimated time: %d seconds', n_combs, time_est);
end

freqinds_low = zeros(n_combs, 2);
freqinds_up = zeros(n_combs, 2);
for i = 1:n_combs
    low = frqs_combs(i,1);
    high = frqs_combs(i,2);
    freqinds_low(i,:) = [find(frqs == low) find(frqs == high - low)]; 
    freqinds_up(i,:) = [find(frqs == low) find(frqs == high)];
end

for proi = 1:nroi 
    for aroi = proi:nroi
        X = data([proi aroi],:,:); % number of regions X epoch length X trails

        % upper freqs
        [BS_up,~] = data2bs_event(X(:,:)', segleng, segshift, epleng, freqinds_up); % BS_up size: 2X2X2
        % normalized by threenorm
        [RTP_up,~] = data2bs_threenorm(X(:,:)', segleng, segshift, epleng, freqinds_up);
        % calculate PAC
        [biv_orig_up, biv_anti_up, biv_orig_up_norm, biv_anti_up_norm] = calc_pac(BS_up, RTP_up);

        % lower freqs
        [BS_low,~] = data2bs_event(X(:,:)', segleng, segshift, epleng, freqinds_low);
        % normalized by threenorm
        [RTP_low,~] = data2bs_threenorm(X(:,:)', segleng, segshift, epleng, freqinds_low);
        % calculate PAC
        [biv_orig_low, biv_anti_low, biv_orig_low_norm, biv_anti_low_norm] = calc_pac(BS_low, RTP_low);

        % PAC_km(f1, f2) = 0.5 * |Bkmm(f1, f2-f1)| + 0.5 * |Bkmm(f1, f2)|
        b_orig(aroi,proi) = mean([biv_orig_up(1) biv_orig_low(1)]); 
        b_orig(proi,aroi) = mean([biv_orig_up(2) biv_orig_low(2)]);
        b_anti(aroi,proi) = mean([biv_anti_up(1) biv_anti_low(1)]);  
        b_anti(proi,aroi) = mean([biv_anti_up(2) biv_anti_low(2)]); 
        
        % normalized versions
        b_orig_norm(aroi,proi) = mean([biv_orig_up_norm(1) biv_orig_low_norm(1)]);
        b_orig_norm(proi,aroi) = mean([biv_orig_up_norm(2) biv_orig_low_norm(2)]);
        b_anti_norm(aroi,proi) = mean([biv_anti_up_norm(1) biv_anti_low_norm(1)]);  
        b_anti_norm(proi,aroi) = mean([biv_anti_up_norm(2) biv_anti_low_norm(2)]);
    end
end