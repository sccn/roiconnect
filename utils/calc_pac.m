function [biv_orig, biv_anti, biv_orig_norm, biv_anti_norm] = calc_pac(BS, RTP)
    % calc_pac - Function to compute bicoherence and antisymmetrized bicoherence
    % Inputs:
    %  BS   - bispectrum frequency data
    %  RTP  - Threenorm data
    % 
    % Outputs:
    %  biv_orig      - Original bicoherence
    %  biv_anti      - Antisymmetrized bicoherence
    %  biv_orig_norm - Normalized original bicoherence
    %  biv_anti_norm - Normalized antisymmetrized bicoherence

    % Calculate bicoherence
    biv_orig = squeeze(([mean(abs(BS(1, 2, 2, :))) mean(abs(BS(2, 1, 1, :)))])); % [Bkmm, Bmkk]
    xx = BS - permute(BS, [2 1 3 4]);
    biv_anti = squeeze(([abs(xx(1, 2, 2, :)) abs(xx(2, 1, 1, :))]));

    % Calculate normalized bicoherence
    bicoh = BS ./ RTP;
    bicoh = mean(bicoh, 4); % average over frequency bands
    biv_orig_norm = ([abs(bicoh(1, 2, 2)) abs(bicoh(2, 1, 1))]);
    xx = bicoh - permute(bicoh, [2 1 3]);
    biv_anti_norm = ([abs(xx(1, 2, 2)) abs(xx(2, 1, 1))]);

end