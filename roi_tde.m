% roi_tde() - Performs time-delay estimation (TDE) based on antisymmetrized bispectra
% between two channels. (cf. Jurhar et al., 2024, Nikias and Pan, 1988).
%
% Usage:
%  EEG = roi_pac(EEG, method, tde_regions, freqbands); 
%
% Inputs:
%   EEG           - EEGLAB dataset with ROI activity computed. Will contain a new EEG.roi.TDE field with the following subfields:
%                     T_XY     - bispectral TDE from region X -> Y
%                     aT_XY    - antisyemmtrized bispectral TDE from region X -> Y
%                     T_YX     - bispectral TDE from region Y -> X
%                     aT_YX    - bispectral TDE from region Y -> X
%                     shift    - time shift (x-axis)
%                     method   - same as 'tde_method'
%                     region_X - region X with its actual label according to EEG.roi.atlas.Scouts
%                     region_Y - region Y with its actual label according to EEG.roi.atlas.Scouts
%   tde_method    - [integer] TDE method, must be between 1:4, open bispectral_TD_est.m for details. Default is 1.
%   tde_regions   - [seed target] Array containing the seed and target region for time-delay estimation. Regions need to be passed as region *indices*, 
%                    e.g., [2 10] will compute time-delays from region 2 -> 10 and 10 -> 2, corresponding to the information flow in both directions separately. 
%                    Default is [] (will throw an error).
%   tde_freqbands - [f1 f2] Array containing the frequency band for which bispectral time-delays will be estimated. Default is [] (broadband).

function EEG = roi_tde(EEG, tde_method, tde_regions, tde_freqbands)

    if nargin < 4
        help roi_tde;
        return
    end

    if ~isfield(EEG, 'roi') || ~isfield(EEG.roi, 'source_roi_data')
        error('Cannot find ROI data - compute ROI data first');
    else
        data = EEG.roi.source_roi_data;
    end

    if isempty(tde_regions)
        help roi_tde
        error('No seed and target regions selected - check the the documentation for the "tde_regions" input parameter in "roi_tde.m"')
    end
    
    % set parameters
    ndat = EEG.pnts * 2;
    segshift = floor(ndat/2);
    segleng = ndat;
    epleng = ndat;
    fres = EEG.srate;
    frqs = sfreqs(2 * fres, EEG.srate);
    maxfreqbins = floor(segleng/2);
    
    % only keeep first PC
    nPCA = EEG.roi.nPCA;
    if nPCA > 1
        warning('Only the first principal component will be used to determine TDE.')
        data = data(1:nPCA:end, :, :);
    end

    % get seed and target region
    X = data(tde_regions(1), :)';
    Y = data(tde_regions(2), :)';

    % compute time-delays
    [T_XY, ~, aT_XY, ~]  = compute_TDE(X, Y, tde_method, tde_freqbands, segleng, segshift, epleng, maxfreqbins, frqs);
    [T_YX, ~, aT_YX, ~]  = compute_TDE(Y, X, tde_method, tde_freqbands, segleng, segshift, epleng, maxfreqbins, frqs);
    
    % store delays in EEG struct
    EEG.roi.TDE.T_XY = T_XY;
    EEG.roi.TDE.aT_XY = aT_XY;
    EEG.roi.TDE.T_YX = T_YX;
    EEG.roi.TDE.aT_YX = aT_YX;
    EEG.roi.TDE.shift = (-floor(segleng/2)+1:floor(segleng/2)-1) / EEG.srate; % x-axis
    EEG.roi.TDE.method = tde_method;
    EEG.roi.TDE.region_X = EEG.roi.atlas.Scouts(tde_regions(1)).Label;
    EEG.roi.TDE.region_Y = EEG.roi.atlas.Scouts(tde_regions(2)).Label;
end

function [T, I, aT, aI] = compute_TDE(X, Y, tde_method, tde_freqbands, segleng, segshift, epleng, maxfreqbins, frqs)

    % compute univariate bispectra
    [B2_xxx] = squeeze(data2bs_univar(X, segleng, segshift, epleng, maxfreqbins));
    para_xyx.chancomb = [1, 2, 1]; 
    [B2_xyx] = data2bs_univar([X, Y], segleng, segshift, epleng, maxfreqbins, para_xyx);
    [B2_yyy] = squeeze(data2bs_univar(Y, segleng, segshift, epleng, maxfreqbins));
    
    % required for antisymmetrization
    para_yxx.chancomb = [2, 1, 1]; 
    [B2_yxx] = data2bs_univar([X, Y], segleng, segshift, epleng, maxfreqbins, para_yxx);

    % compute time-delays
    if isempty(tde_freqbands)
        [T, I] = bispectral_TD_est(B2_xxx, B2_xyx, B2_yyy, tde_method, [], 1);
        [aT, aI] = bispectral_TD_est(B2_xxx, B2_xyx - B2_yxx, B2_yyy, tde_method, [], 1); % antisymmetrization
    else
        fmask = true(1, length(frqs));
        fmask(frqs < tde_freqbands(1) | frqs > tde_freqbands(2)) = 0;
        [T, I] = bispectral_TD_est(B2_xxx, B2_xyx, B2_yyy, tde_method, fmask(1:end-1), 1);
        [aT, aI] = bispectral_TD_est(B2_xxx, B2_xyx - B2_yxx, B2_yyy, tde_method, fmask(1:end-1), 1); % antisymmetrization
    end
end