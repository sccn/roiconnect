% (C) Stefan Haufe 2018

%% settings

datasets = {'NDARAA948VFH'};

for isb = 1:length(datasets)
    
    dataset = datasets{isb};
    
    inputdir = ['analysis_output/sub_', dataset, '/preprocessing/data/'];
    
    %%% Creating result folder
    resultsdir = ['analysis_output/sub_', dataset, '/source_connectivity/'];
    mkdir([resultsdir 'data/'])
    
    brainstormprojectdir = 'Brainstorm/ChildMindICBM152/';
    
    % in this script, a standard head model (ICBM152), computed as a 3-shell
    % BEM with OpenMEEG in Brainstorm is used. This can be replaced with an
    % individual head model computed in Brainstorm by just exchanging the
    % respective paths to the BS database.
    
    % load the cortical surface for which the leadfield has been computed. I
    % used a downsampled version with only 2000 nodes, and I used a surface
    % My idea is to use an ROI based analyses, so the intial resolution should
    % not be crucial. But one could go up to about 5000 nodes.
    cortex = load([brainstormprojectdir 'anat/@default_subject/tess_cortex_mid_low_2000V.mat']);
    nvox = size(cortex.Vertices, 1);
    
    % load a very high resolution version of the same surface for plotting
    cortex_highres = load([brainstormprojectdir 'anat/@default_subject/tess_cortex_mid_high.mat']);
    
    % make brainstorm coordinate system consistent with MNI coordinates for
    % plotting (in terms of axis directions)
    cortex.Vertices = cortex.Vertices(:, [2 1 3]);
    % in particular, the left-right axis needs to be flipped
    cortex.Vertices(:, 1) = -cortex.Vertices(:, 1);
    
    % same for highres_cortex
    cortex_highres.Vertices = cortex_highres.Vertices(:, [2 1 3]);
    cortex_highres.Vertices(:, 1) = -cortex_highres.Vertices(:, 1);
    
    % calculate extrapolation from coarse to highres cortex, takes one minute
    mi = []; in_low_to_high = [];
    for ii = 1:size(cortex_highres.Vertices, 1);
        %   ii
        [mi(ii) in_low_to_high(ii)] = min(eucl(cortex_highres.Vertices(ii, :), cortex.Vertices));
    end
    
    % load a 3-shell BEM model for the ICBM152 head with EGI Hydrocel 129 cap
    % that I precomputed manually in Brainstorm (I can show you the necessary
    % steps if needed)
    headmodel = load([brainstormprojectdir 'data/@default_study/headmodel_surf_openmeeg.mat']);
    
    % make format compatible with my routines
    leadfield = permute(reshape(headmodel.Gain, [], 3, nvox), [1 3 2]);
    
    
    % model order for the spectral/connectivity analysis
    morder = 20;
    
    % regularization parameter for eloreta
    eloreta_reg = 0.05;
    
    % how many PCA components to keep for each ROI
    nPCA = 3;
    
    % colormap
    load cm17;
    
    % plot smooth cortical surfaces
    smooth_cortex = 0.35;
    
    EEG = pop_loadset([inputdir '/prep_interp.set']);
    
    % use frequency resolution of 0.5 Hz
    fres = EEG.srate;
    
    % from the MVGC toolbox, compute frequencies in Hz for a
    frqs = sfreqs(fres, EEG.srate);
    
    % frequency bin indices for some EEG rhythms
    bands = {[4 8], [8 13], [13 28]};
    bandnames = {'theta', 'alpha', 'beta'};
    for iba = 1:length(bands)
        frq_inds{iba} = find(frqs >= bands{iba}(1) & frqs < bands{iba}(2));
    end
    
    % number of bootstrap samples
    nbootstrap = 1;
    
    %% source reconstruction
    
    % common average reference transform
    H = eye(EEG.nbchan) - ones(EEG.nbchan) ./ EEG.nbchan;
    
    % apply to data and leadfield
    EEG.data = reshape(H*EEG.data(:, :), EEG.nbchan, EEG.pnts, EEG.trials);
    leadfield = reshape(H*leadfield(:, :), EEG.nbchan, nvox, 3);
    
    % eLORETA inverse projection kernel
    P_eloreta = mkfilt_eloreta_v2(leadfield, eloreta_reg);
    
    % project to source space
    source_voxel_data = reshape(EEG.data(:, :)'*P_eloreta(:, :), EEG.pnts*EEG.trials, nvox, 3);
    
    % plot broadband source power per voxel
    P = reshape(sum(sum(source_voxel_data.^2, 1), 3), [], 1);
    P_dB = 10*log10(P);
%     allplots_cortex_BS(cortex_highres, P_dB(in_low_to_high), [min(P_dB) max(P_dB)], cm17a, 'power [dB]', smooth_cortex, ...
%         [resultsdir 'print/eloreta' num2str(eloreta_reg) '_voxel_broadband_power']);
    allplots_cortex_BS(cortex, P_dB, [min(P_dB) max(P_dB)], cm17a, 'power [dB]', smooth_cortex, ...
        [resultsdir 'print/eloreta' num2str(eloreta_reg) '_voxel_broadband_power']);
    
    % % or plot the same just on the lowres cortex
    % allplots_cortex_BS(cortex, P, [min(P) max(P)], cm17a, 'power [a.u.]', smooth_cortex, ...
    %   ['print/eloreta' num2str(eloreta_reg) '_voxel_broadbandPower']);
    
    % number of ROIs in the Desikan-Killiany Atlas
    nROI = length(cortex.Atlas(3).Scouts);
    
    % ROI labels
    labels = {cortex.Atlas(3).Scouts.Label};
    
    % keep only the first nPCA strongest components for each ROI
    source_roi_data = [];
    for iROI = 1:nROI
        ind_roi = cortex.Atlas(3).Scouts(iROI).Vertices;
        data_  = source_voxel_data(:, ind_roi, :);
        source_roi_power(iROI) = sum(var(data_(:, :)))';
        source_roi_power_norm(iROI) = source_roi_power(iROI)/length(ind_roi);
        
        % optional z-scoring, this makes the PCA independent of the power in each
        % voxel, and favors to find components that are consistently expressed in
        % many voxels rather than only in a few voxels with strong power (which
        % may leak from a neighboring region)
        data_(:, :) = zscore(data_(:, :));
        
        [data_, ~, ~] = svd(data_(:, :), 'econ');
        source_roi_data(:, :, iROI) = data_(:, 1:nPCA);
    end
    
    % version with nPCA components
    source_roi_data = permute(reshape(source_roi_data, EEG.pnts, EEG.trials, []), [3 1 2]);
    
    % plot broadband source power per region
    source_roi_power_norm_dB = 10*log10(source_roi_power_norm);
    allplots_cortex_BS(cortex_highres, source_roi_power_norm_dB, [min(source_roi_power_norm_dB) max(source_roi_power_norm_dB)], ...
        cm17a, 'power [dB]', smooth_cortex, ...
        [resultsdir 'print/eloreta' num2str(eloreta_reg) '_roi_broadband_power']);
    
    %% spectral and connectivity analysis
    
    % to test TRGC between ROIs (that is, pairs of nPCA-dimensional spaces), we
    % need to compute these indices
    inds = {}; ninds = 0;
    for iroi = 1:nROI
        for jroi = (iroi+1):nROI
            inds{ninds+1} = {(iroi-1)*nPCA + [1:nPCA], (jroi-1)*nPCA + [1:nPCA]};
            inds{ninds+2} = {(jroi-1)*nPCA + [1:nPCA], (iroi-1)*nPCA + [1:nPCA]};
            ninds = ninds + 2;
        end
    end
    
    % compute time reversed spectral Granger causality between all pairs of ROIs
    TRGC = data2sctrgc(source_roi_data, fres, morder, 0, nbootstrap, [], inds);
    
    % calculation of net TRGC scores (i->j minus j->i), recommended procedure
    TRGCnet = TRGC(:, 1:2:end)-TRGC(:, 2:2:end);
    
    % create a ROI x ROI connectivity matrix, if needed
    % TRGCmat(f, ii, jj) is net TRGC from jj to ii
    TRGCmat = [];
    iinds = 0;
    for iroi = 1:nROI
        for jroi = (iroi+1):nROI
            iinds = iinds + 1;
            TRGCmat(:, iroi, jroi) = -TRGCnet(:, iinds);
            TRGCmat(:, jroi, iroi) = TRGCnet(:, iinds);
        end
    end
    
    
    % plot Roi to ROI TRGC across entire spectrum
    figure; imagesc(squeeze(mean(TRGCmat))); colorbar
    
    % the same just for the alpha band
    figure; imagesc(squeeze(mean(TRGCmat(frq_inds{2}, :, :)))); colorbar
    
    % plot global net TRGC in alpha band on cortex,
    % "global net" is defined here as the sum over the interactions to all other brain
    % areas
    atrgc = squeeze(mean(mean(TRGCmat(frq_inds{2}, :, :), 1), 2));
    
    % red = net sender, blue = net receiver
    allplots_cortex_BS(cortex_highres, atrgc, [-max(abs(atrgc)) max(abs(atrgc))], cm17, 'TRGC', smooth_cortex, ...
        [resultsdir 'print/eloreta' num2str(eloreta_reg) '_roi_alpha_TRGCglobalnet']);
    
    
    % compute cross-spectrum, takes a while
    conn = data2spwctrgc(source_roi_data, fres, morder, 0, nbootstrap, [], {'CS'});
    
    % power spectral density from the cross-spectrum
    PS = cs2psd(conn.CS);
    
    figure; semilogy(frqs, PS); grid on
    
    % plot alpha power on cortex
    apow = squeeze(sum(sum(reshape(PS(frq_inds{2}, :), [], nPCA, nROI), 1), 2)).*source_roi_power_norm';
    
    apow_dB = 10*log10(apow);
    allplots_cortex_BS(cortex_highres, apow_dB, [min(apow_dB) max(apow_dB)], cm17a, 'power [dB]', smooth_cortex, ...
        [resultsdir 'print/eloreta' num2str(eloreta_reg) '_roi_alpha_power']);
    
    % absolute value of the imaginary part of coherence, a robust
    % non-directional measure of synchronization with non-zero phase lag
    absiCOH = abs(imag(cs2coh(conn.CS)));
    absiCOH = squeeze(mean(mean(reshape(absiCOH, fres+1, 3, nROI, 3, nROI), 2), 4));
    
    figure; imagesc(squeeze(mean(absiCOH(frq_inds{2}, :, :)))); colorbar
    
    % global net absolute value of the imaginary part of coherence =
    % how much interaction each ROI is involved in
    netabsiCOH = mean(squeeze(mean(absiCOH(frq_inds{2}, :, :))), 2);
    
    % red = highly interacting ROI
    allplots_cortex_BS(cortex_highres, netabsiCOH, [min(netabsiCOH) max(netabsiCOH)], cm17a, 'net |iCOH|', smooth_cortex, ...
        [resultsdir 'print/eloreta' num2str(eloreta_reg) '_roi_alpha_abs_iCOH_globalnet']);
    
    save([resultsdir 'data/TRGC'], 'TRGC', 'conn')
    
end

