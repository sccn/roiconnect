% pac_bispec() - Estimates phase-amplitude coupling (PAC) by computing  
%                bicoherence. Direct relations between PAC and bispectra 
%                are shown in Zandvoort and Nolte, 2021.
%
% Inputs:
%  data   - nchan x ntimepoints x ntrials ROI data, requires prior source reconstruction 
%  params - Set containing parameters for bispectrum computation
%  
% Output:
%   b_orig - ROI x ROI bispectrum
%   b_anti - ROI x ROI antisymmetrized bispectrum
%   b_orig_norm - ROI x ROI normalized bispectrum
%   b_anti_norm - ROI x ROI normalized antisymmetrized bispectrum
%
% Copyright (C) Franziska Pellegrini, franziska.pellegrini@charite.de

function [b_orig, b_anti, b_orig_norm,b_anti_norm] = pac_bispec(data, params)

nroi = params.nROI;
segleng = params.segleng;
segshift = params.segshift;
epleng = params.epleng;
filt = params.filt;
fs = params.fs;

fres = fs;
frqs = sfreqs(fres, fs);
freqinds_low = [find(frqs==filt.low) find(frqs==filt.high-filt.low)];
freqinds_up = [find(frqs==filt.low) find(frqs==filt.high)];

for proi = 1:nroi
    for aroi = proi:nroi
        
        % upper freqs
        X = data([proi aroi],:,:); 
        [bs_up,~] = data2bs_event(X(:,:)', segleng, segshift, epleng, freqinds_up);
        biv_orig_up = ([abs(bs_up(1, 2, 2)) abs(bs_up(2, 1, 1))]);
        xx = bs_up - permute(bs_up, [2 1 3]); %Bkmm - Bmkm
        biv_anti_up = ([abs(xx(1, 2, 2)) abs(xx(2, 1, 1))]);
        
        % normalized by threenorm
        [RTP_up,~]=data2bs_threenorm(X(:,:)',segleng,segshift,epleng,freqinds_up);
        bicoh_up=bs_up./sqrt(RTP_up);
        biv_orig_up_norm = ([abs(bicoh_up(1, 2, 2)) abs(bicoh_up(2, 1, 1))]);
        xx=bicoh_up-permute(bicoh_up, [2 1 3]);
        biv_anti_up_norm = ([abs(xx(1, 2, 2)) abs(xx(2, 1, 1))]);
        
        % lower freqs
        [bs_low,~] = data2bs_event(X(:,:)', segleng, segshift, epleng, freqinds_low);
        biv_orig_low = ([abs(bs_low(1, 2, 2)) abs(bs_low(2, 1, 1))]);
        xx = bs_low - permute(bs_low, [2 1 3]);
        biv_anti_low =([abs(xx(1, 2, 2)) abs(xx(2, 1, 1))]);
        
        % normalized by threenorm
        [RTP_low,~]=data2bs_threenorm(X(:,:)',segleng,segshift,epleng,freqinds_low);
        bicoh_low=bs_low./sqrt(RTP_low);
        biv_orig_low_norm = ([abs(bicoh_low(1, 2, 2)) abs(bicoh_low(2, 1, 1))]);
        xx=bicoh_low-permute(bicoh_low, [2 1 3]);
        biv_anti_low_norm = ([abs(xx(1, 2, 2)) abs(xx(2, 1, 1))]);
        
        b_orig(aroi,proi) = mean([biv_orig_up(1) biv_orig_low(1)]); % mean across the two possible bispec combinations 
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