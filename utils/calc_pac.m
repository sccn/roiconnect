% calc_pac - Function to compute absolute values of the cross-bispectrum/bicoherence and the antisymmetrized 
%            cross-bispectrum/bicoherence, such that biv_orig corresponds to |B(f1,f2)|. Used for the final 
%            computation of PAC_km(f1, f2) = 0.5 * |Bkmm(f1, f2-f1)| + 0.5 * |Bkmm(f1, f2)|.
%
% Inputs:
%  BS   - [array] (nchan x nchan x nchan x nfreqpairs x nshuf) Cross-bispectral tensor
%  RTP  - [array] (nchan x nchan x nchan) Threenorm tensor to normalize the cross-bispectrum, used to compute cross-bicoherence
% 
% Outputs:
%  biv_orig      - [array] (nchan x nshuf) Absolute value of the original cross-bispectrum, corresponding to [Bkmm(f1, f2), Bmkk(f1, f2)]
%  biv_anti      - [array] (nchan x nshuf) Absolute value of the antisymmetrized cross-bispectrum
%  biv_orig_norm - [array] (nchan x nshuf) Absolute value of the normalized original cross-bispectrum (cross-bicoherence)
%  biv_anti_norm - [array] (nchan x nshuf) Absolute value of the normalized antisymmetrized cross-bispectrum (antisymmetrized cross-bicoherence)
%
% Authors:
%   Tien Dung Nguyen, tien-dung.nguyen@charite.de,  
%   Zixuan Liu, zixuan.liu@campus.tu-berlin.de

function [biv_orig, biv_anti, biv_orig_norm, biv_anti_norm] = calc_pac(BS, RTP)              

    % Calculate absolute values of the cross-bispectrum
    biv_orig = squeeze(([abs(mean(BS(1, 2, 2, :, :), 4)) abs(mean(BS(2, 1, 1, :, :), 4))])); % [Bkmm, Bmkk], average over frequency bands (4th dimension)
    xx = BS - permute(BS, [2 1 3 4 5]); 
    biv_anti = squeeze(([abs(mean(xx(1, 2, 2, :,:), 4)) abs(mean(xx(2, 1, 1, :, :), 4))]));

    % Calculate absolute values of the cross-bicoherence
    bicoh = BS ./ RTP;
    biv_orig_norm = squeeze(([abs(mean(bicoh(1, 2, 2, :, :), 4)) abs(mean(bicoh(2, 1, 1, :, :), 4))])); % [Bkmm, Bmkk], average over frequency bands
    xx = bicoh - permute(bicoh, [2 1 3 4 5]);
    biv_anti_norm = squeeze(([abs(mean(xx(1, 2, 2, :, :), 4)) abs(mean(xx(2, 1, 1,:, :), 4))]));

end