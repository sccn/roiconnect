function [biv_orig, biv_anti, biv_orig_norm, biv_anti_norm] = calc_pac(BS, RTP)
    % calc_pac - Function to compute bicoherence and antisymmetrized bicoherence
    % Inputs:
    %  BS   - bispectrum frequency data
    %         nchan by nchan by nchan by number_of_frequency_pairs (by number_of_segments) tensor such that 
    %         BS(i,j,k,f)=<x(f1)_i*x(f2)_j*conj(x(f1+f2-1)_k)>
    %         where f1=freqpairs(f,1) and  f2=freqpairs(f,2), x=fft(data) and the average is over epeochs and segments

    %  RTP  - Threenorm data
    %         nchan by nchan by nchan by number_of_frequency_pairs
    %         (by number_of_segments) tensor
    % 
    % Outputs:
    %  biv_orig      - Original bicoherence
    %  biv_anti      - Antisymmetrized bicoherence
    %  biv_orig_norm - Normalized original bicoherence
    %  biv_anti_norm - Normalized antisymmetrized bicoherence
    %
    % Authors:
    %   Zixuan Liu, zixuan.liu@campus.tu-berlin.de
    %   Tien Dung Nguyen, tien-dung.nguyen@charite.de                

    % Calculate bicoherence
    biv_orig = squeeze(([mean(abs(BS(1, 2, 2, :, :))) mean(abs(BS(2, 1, 1, :, :)))])); % [Bkmm, Bmkk] % add dimension
    xx = BS - permute(BS, [2 1 3 4 5]);
    biv_anti = squeeze(([abs(xx(1, 2, 2, :,:)) abs(xx(2, 1, 1, :,:))]));

    % Calculate normalized bicoherence
    bicoh = BS ./ RTP;
    bicoh = mean(bicoh, 4); % average over frequency bands
    biv_orig_norm = ([abs(bicoh(1, 2, 2)) abs(bicoh(2, 1, 1))]);
    xx = bicoh - permute(bicoh, [2 1 3 4 5]);
    biv_anti_norm = ([abs(xx(1, 2, 2, :, :)) abs(xx(2, 1, 1,:, :))]);

end