function [conn, nlags] = data2spwctrgc(data, fres, nlags, cond, nboot, maxfreq, output, verbose, varargin)
% Epoched time series data to spectral pairwise-conditional time-reversed Granger causality
% and other useful quantities (cross-spectrum, (time reversed) directed
% transfer function, (time reversed) partial directed cCOHerence)
%
% conn.measures(:, ii, jj, :) refers to the flow from jj->ii
%
% (C) 2018 Stefan Haufe
%
% This code implements/calls methods presented in the follow papers. Please
% cite the appropriate ones if you use it for a publication.
%
% Haufe, S., Nikulin, V. V., M?ller, K. R., & Nolte, G. (2013). A critical
% assessment of connectivity measures for EEG data: a simulation study.
% Neuroimage, 64, 120-133.
%
% Haufe, S., Nikulin, V. V., & Nolte, G. (2012, March). Alleviating the
% influence of weak data asymmetries on granger-causal analyses.
% In International Conference on Latent Variable Analysis and Signal
% Separation (pp. 25-33). Springer, Berlin, Heidelberg.
%
% Winkler, I., Panknin, D., Bartz, D., M?ller, K. R., & Haufe, S. (2016).
% Validity of time reversal for testing Granger causality. IEEE Transactions
% on Signal Processing, 64(11), 2746-2760.
%
% L Faes, D Marinazzo, S Stramaglia, 'Multiscale information decomposition:
% exact computation for multivariate Gaussian processes', Entropy, special
% issue on Multivariate entropy measures and their applications, 2017, 19(8),
% 408. DOI: 10.3390/e19080408.
%
% Barnett, L.; Seth, A.K. Granger causality for state-space models.
% Phys. Rev. E 2015, 91, 040101.
%
% Solo, V. State-space analysis of Granger-Geweke causality measures with
% application to fMRI. Neural Computation 2016, 28, 914?949.
%
% Barnett, L., & Seth, A. K. (2014). The MVGC multivariate Granger causality
% toolbox: a new approach to Granger-causal inference. Journal of
% neuroscience methods, 223, 50-68.

% decode input parameters
% -----------------------
g = finputcheck(varargin, { ...
    'freqresolution'  'integer'  { }    0;
    'roi_selection'   'cell'     { }    { }; ...
    'nPCA'            'integer'  { }    3  }, 'data2spwctrgc'); 
if ischar(g), error(g); end

% choose ROIs if desired, take number of PCs into account
if ~isempty(g.roi_selection)
    nROI = length(g.roi_selection);
    nPCA = g.nPCA;
    data_new = zeros(nROI * nPCA, size(data, 2), size(data, 3));
    
    start_idx_new  = 1;
    end_idx_new = nPCA;
    for iroi = 1:nROI
        end_idx = g.roi_selection{iroi} * nPCA;
        start_idx = end_idx - (nPCA - 1);
        data_new(start_idx_new:end_idx_new, :, :) = data(start_idx:end_idx, :, :);

        start_idx_new = start_idx_new + nPCA;
        end_idx_new = start_idx_new + nPCA - 1;
    end
    data = data_new;
end

[nchan, ndat, nepo] = size(data);

ninds = nchan*(nchan-1);

if nargin < 2 || isempty(fres)
    fres = ndat;
end
freqs = linspace(0, 1, fres+1);

if nargin < 3 || isempty(nlags)
    nlags = -10;
end

if nargin < 4 || isempty(cond)
    cond = 0;
end

if nargin < 5 || isempty(nboot)
    nboot = 0;
else
    if nepo < 2
        warning('Not enough epochs to perform bootstrap.')
    end
end

if nargin < 6 || isempty(maxfreq)
    maxfreq = fres+1;
end

if nargin < 7 || isempty(output)
    output = {'TRGC'};
end

if nargin < 8 || isempty(verbose)
    verbose = 0;
end

freqs = freqs(1:maxfreq);
z = exp(-1i*pi*freqs);

if  nlags < 0
    if cond
        [~,~,~, nlags] = tsdata_to_infocrit(data, -nlags, [], 0);
    else
        nlagsall = [];
        for ii = 1:min(nchan*(nchan-1)/2, -nlags)
            pe = randperm(nchan);
            [~,~,~, nlags_] = tsdata_to_infocrit(data(pe(1:2), :, :), -nlags, [], 0);
            nlagsall = [nlagsall nlags_];
        end
        nlags = floor(median(nlagsall));
    end
end

CSpara = [];
CSpara.subave = 0;
CSpara.mywindow = hanning(ndat)./sqrt(hanning(ndat)'*hanning(ndat));
CSpara.freqresolution = g.freqresolution;

clear TRGC GC TRPDC PDC TRDTF DTF CS aCOH iCOH cCOH wPLI

if abs(nboot) < 1 % no bootstrap

    %     % data to autocovariance
    %     G = tsdata_to_autocov(data, nlags);

    if verbose; disp(['computing cross-spectrum']); end

    [CS, ~, wPLI, ~] = data2cs_event(data(:, :)', ndat, floor(ndat/2), ndat, [], CSpara);
    ASR = log(norm(imag(CS(:)))/norm(real(CS(:))));
    fprintf('ASR: %1.4f\n', ASR);

    maxfreq = size(CS,3);

    if ~isempty(intersect(output, {'cCOH', 'iCOH', 'aCOH'}))
        clear cCOH iCOH aCOH
        for ifreq = 1:maxfreq
            clear pow
            pow = real(diag(CS(:,:,ifreq)));
            cCOH(:,:,ifreq) = CS(:,:,ifreq) ./ sqrt(pow*pow');
            aCOH(:,:,ifreq) = abs(CS(:,:,ifreq)) ./ sqrt(pow*pow');
            iCOH(:,:,ifreq) = abs(imag(CS(:,:,ifreq) ./ sqrt(pow*pow')));
        end
    end

    if ~isempty(intersect(output, {'TRGC', 'GC', 'TRPDC', 'PDC', 'TRDTF', 'DTF'}))
        G = cpsd_to_autocov(CS, nlags);
    end

    if isempty(intersect(output, {'CS'}))
        clear CS
    end

    if ~isempty(intersect(output, {'TRGC', 'GC', 'TRPDC', 'PDC', 'TRDTF', 'DTF'}))

        if cond
            % (time-reversed) GC conditioned on all other variables

            % autocovariance to full forward VAR model
            [A, SIG] = autocov_to_var3(G);

            % forward VAR model to state space VARMA models
            [eA2, eC2, eK2, eV2, eVy2] = varma2iss(reshape(A, nchan, []), [], SIG, eye(nchan));

            if ~isempty(intersect(output, {'TRGC', 'TRPDC', 'TRDTF'}))
                % backward autocovariance to full backward VAR model
                [AR, SIGR] = autocov_to_var3(permute(G, [2 1 3]));

                % backward VAR to VARMA
                [eA2R, eC2R, eK2R, eV2R, eVy2R] = varma2iss(reshape(AR, nchan, []), [], SIGR, eye(nchan));
            end

            if ~isempty(intersect(output, {'TRGC', 'GC'}))
                %% loop over sender/receiver combinations to compute (time-reversed) GC
                iind = 0;
                fprintf('Progress of %d:', nchan)
                for ii = 1:nchan
                    if mod(ii,10) == 0
                        fprintf('%d', ii);
                    else
                        fprintf('.');
                    end
                    for jj = (ii+1):nchan
                        iind = iind + 1;
                        if verbose; disp(['testing connection ' num2str(iind) '/' num2str(ninds) ': [' num2str(ii) '] -> [' num2str(jj) '], conditional']); end

                        GC(:, ii, jj) = iss_SGC(eA2, eC2, eK2, eV2, z, ii, jj);

                        if ~isempty(intersect(output, {'TRGC'}))
                            GCR = iss_SGC(eA2R, eC2R, eK2R, eV2R, z, ii, jj);
                            TRGC(:, ii, jj) = GC(:, ii, jj) - GCR';
                        end

                        iind = iind + 1;
                        if verbose; disp(['testing connection ' num2str(iind) '/' num2str(ninds) ': [' num2str(jj) '] -> [' num2str(ii) '], conditional']); end

                        GC(:, jj, ii) = iss_SGC(eA2, eC2, eK2, eV2, z, jj, ii);

                        if ~isempty(intersect(output, {'TRGC'}))
                            GCR = iss_SGC(eA2R, eC2R, eK2R, eV2R, z, jj, ii);
                            TRGC(:, jj, ii) = GC(:, jj, ii) - GCR';
                        end
                    end
                end
                fprintf('\n');
            end

            % PDC
            if ~isempty(intersect(output, {'TRPDC', 'PDC'}))
                PDC = permute(abs(est_PDC(A, freqs/2)), [3 1 2]);
            end

            % TRPDC
            if ~isempty(intersect(output, {'TRPDC'}))
                TRPDC = PDC - permute(abs(est_PDC(AR, freqs/2)), [3 1 2]);
            end

            % DTF
            if ~isempty(intersect(output, {'TRDTF', 'DTF'}))
                DTF = permute(abs(est_DTF(A, SIG, freqs/2)), [3 1 2]);
            end

            % TRDTF
            if ~isempty(intersect(output, {'TRDTF'}))
                TRDTF = DTF - permute(abs(est_DTF(AR, SIGR, freqs/2)), [3 1 2]);
            end

        else
            % (time-reversed) GC just between sender and receiver

            % loop over sender/receiver combinations to compute time-reversed GC
            iind = 0;
            fprintf('Progress of %d:', nchan)
            for ii = 1:nchan
                if mod(ii,10) == 0
                    fprintf('%d', ii);
                else
                    fprintf('.');
                end
                for jj = (ii+1):nchan
                    iind = iind + 1;
                    if verbose; disp(['testing connection ' num2str(iind) '/' num2str(ninds) ': [' num2str(ii) '] -> [' num2str(jj) ']']); end

                    subset = [ii jj];

                    % autocovariance to full forward VAR model
                    [A, SIG] = autocov_to_var3(G(subset, subset, :));

                    % forward VAR model to state space VARMA models
                    [eA2, eC2, eK2, eV2, eVy2] = varma2iss(reshape(A, 2, []), [], SIG, eye(2));

                    if ~isempty(intersect(output, {'TRGC', 'TRPDC', 'TRDTF'}))
                        % backward autocovariance to full backward VAR model
                        [AR, SIGR] = autocov_to_var3(permute(G(subset, subset, :), [2 1 3]));

                        % backward VAR to VARMA
                        [eA2R, eC2R, eK2R, eV2R, eVy2R] = varma2iss(reshape(AR, 2, []), [], SIGR, eye(2));
                    end

                    if ~isempty(intersect(output, {'TRGC', 'GC'}))
                        % GC and TRGC computation
                        GC(:, ii, jj) = iss_SGC(eA2, eC2, eK2, eV2, z, 1, 2);

                        if ~isempty(intersect(output, {'TRGC'}))
                            GCR = iss_SGC(eA2R, eC2R, eK2R, eV2R, z, 1, 2);
                            TRGC(:, ii, jj) = GC(:, ii, jj) - GCR';
                        end


                        iind = iind + 1;
                        if verbose; disp(['testing connection ' num2str(iind) '/' num2str(ninds) ': [' num2str(jj) '] -> [' num2str(ii) ']']); end

                        % GC and TRGC computation
                        GC(:, jj, ii) = iss_SGC(eA2, eC2, eK2, eV2, z, 2, 1);

                        if ~isempty(intersect(output, {'TRGC'}))
                            GCR = iss_SGC(eA2R, eC2R, eK2R, eV2R, z, 2, 1);
                            TRGC(:, jj, ii) = GC(:, jj, ii) - GCR';
                        end
                    end

                    % PDC
                    if ~isempty(intersect(output, {'TRPDC', 'PDC'}))
                        PDC_ = abs(est_PDC(A, freqs/2));
                        PDC(:, ii, jj) = PDC_(1, 2, :);
                        PDC(:, jj, ii) = PDC_(2, 1, :);
                    end

                    % TRPDC
                    if ~isempty(intersect(output, {'TRPDC'}))
                        TRPDC_ = PDC_ - abs(est_PDC(AR, freqs/2));
                        TRPDC(:, ii, jj) = TRPDC_(1, 2, :);
                        TRPDC(:, jj, ii) = TRPDC_(2, 1, :);
                    end

                    % DTF
                    if ~isempty(intersect(output, {'TRDTF', 'DTF'}))
                        DTF_ = abs(est_DTF(A, SIG, freqs/2));
                        DTF(:, ii, jj) = DTF_(1, 2, :);
                        DTF(:, jj, ii) = DTF_(2, 1, :);
                    end

                    % TRDTF
                    if ~isempty(intersect(output, {'TRDTF'}))
                        TRDTF_ = DTF_ - abs(est_DTF(AR, SIGR, freqs/2));
                        TRDTF(:, ii, jj) = TRDTF_(1, 2, :);
                        TRDTF(:, jj, ii) = TRDTF_(2, 1, :);
                    end

                end
            end
            fprintf('\n');
        end
    end

else % bootstrap

    % loop over bootstrap samples
    for iboot = 1:nboot
        bootinds = randi(nepo, nepo, 1);

        %       % data to autocovariance
        %       G = tsdata_to_autocov(data(:, :, bootinds), nlags);

        if verbose; disp(['bootstrap run ' num2str(iboot) '/' num2str(nboot) ', computing cross-spectrum']); end

        if ~isempty(intersect(output, {'wPLI'}))
            [CS_, wPLI_] = tsdata_to_cpsdwpli(data(:, :, bootinds), fres, 'MT', ndat, [], 1);

            wPLI_(isnan(wPLI_)) = 0;
            wPLI_ = permute(wPLI_, [3, 1, 2]);
            wPLI(:, :, :, iboot) = wPLI_(1:maxfreq, :, :);
        else
            %         CS_ = tsdata_to_cpsd(data(:, :, bootinds), fres, 'MT', ndat, [], 1);
            data_ = data(:, :, bootinds);
            CS_ = data2cs_event(data_(:, :)', ndat, ndat, ndat, fres+1, CSpara);
        end

        if ~isempty(intersect(output, {'CS'}))
            CS(:, :, :, iboot) = CS_(:, :, 1:maxfreq);
        end

        if ~isempty(intersect(output, {'cCOH', 'aCOH', 'iCOH'}))
            clear cCOH_ aCOH_ iCOH_
            for ifreq = 1:maxfreq
                clear pow
                pow = real(diag(CS_(:,:,ifreq)));
                cCOH_(:,:,ifreq) = CS_(:,:,ifreq) ./ sqrt(pow*pow');
                aCOH_(:,:,ifreq) = abs(CS_(:,:,ifreq)) ./ sqrt(pow*pow');
                iCOH_(:,:,ifreq) = abs(imag(CS_(:,:,ifreq) ./ sqrt(pow*pow')));
            end
            if ~isempty(intersect(output, {'cCOH', 'aCOH', 'COH'}))
                cCOH(:, :, :, iboot) = cCOH_;
                aCOH(:, :, :, iboot) = aCOH_;
                iCOH(:, :, :, iboot) = iCOH_;
            end
        end

        if ~isempty(intersect(output, {'GC', 'TRGC'}))
            G = cpsd_to_autocov(CS_, nlags);
        end

        clear CS_

        if ~isempty(intersect(output, {'TRGC', 'GC', 'TRPDC', 'PDC', 'TRDTF', 'DTF'}))

            if cond
                % (time-reversed) GC conditioned on all other variables

                % autocovariance to full forward VAR model
                [A, SIG] = autocov_to_var3(G);

                % forward VAR model to state space VARMA models
                [eA, eC, eK, eV, eVy] = varma2iss(reshape(A, nchan, []), [], SIG, eye(nchan));

                if ~isempty(intersect(output, {'TRGC', 'TRPDC', 'TRDTF'}))
                    % backward autocovariance to full backward VAR model
                    [AR, SIGR] = autocov_to_var3(permute(G, [2 1 3]));

                    % backward VAR to VARMA
                    [eAR, eCR, eKR, eVR, eVyR] = varma2iss(reshape(AR, nchan, []), [], SIGR, eye(nchan));
                end

                if ~isempty(intersect(output, {'TRGC', 'GC'}))
                    % loop over sender/receiver combinations to compute (time-reversed) GC
                    iind = 0;
                    for ii = 1:nchan
                        for jj = (ii+1):nchan
                            iind = iind + 1;
                            if verbose; disp(['bootstrap run ' num2str(iboot) '/' num2str(nboot) ', testing connection ' num2str(iind) '/' num2str(ninds) ': [' num2str(ii) '] -> [' num2str(jj) '], conditional']); end

                            GC(:, ii, jj, iboot) = iss_SGC(eA, eC, eK, eV, z, ii, jj);

                            if ~isempty(intersect(output, {'TRGC'}))
                                GCR = iss_SGC(eAR, eCR, eKR, eVR, z, ii, jj);
                                TRGC(:, ii, jj, iboot) = GC(:, ii, jj, iboot) - GCR';
                            end

                            iind = iind + 1;
                            if verbose; disp(['bootstrap run ' num2str(iboot) '/' num2str(nboot) ', testing connection ' num2str(iind) '/' num2str(ninds) ': [' num2str(jj) '] -> [' num2str(ii) '], conditional']); end

                            GC(:, jj, ii, iboot) = iss_SGC(eA, eC, eK, eV, z, jj, ii);

                            if ~isempty(intersect(output, {'TRGC'}))
                                GCR = iss_SGC(eAR, eCR, eKR, eVR, z, jj, ii);
                                TRGC(:, jj, ii, iboot) = GC(:, jj, ii, iboot) - GCR';
                            end

                        end
                    end
                end

                % PDC
                if ~isempty(intersect(output, {'TRPDC', 'PDC'}))
                    PDC(:, :, :, iboot) = permute(abs(est_PDC(A, freqs/2)), [3 1 2]);
                end

                % TRPDC
                if ~isempty(intersect(output, {'TRPDC', 'PDC'}))
                    TRPDC(:, :, :, iboot) = PDC(:, :, :, iboot) - permute(abs(est_PDC(AR, freqs/2)), [3 1 2]);
                end

                % DTF
                if ~isempty(intersect(output, {'TRDTF', 'DTF'}))
                    DTF(:, :, :, iboot) = permute(abs(est_DTF(A, SIG, freqs/2)), [3 1 2]);
                end

                % TRDTF
                if ~isempty(intersect(output, {'TRDTF'}))
                    TRDTF(:, :, :, iboot) = DTF(:, :, :, iboot) - permute(abs(est_DTF(AR, SIGR, freqs/2)), [3 1 2]);
                end

            else
                % (time-reversed) GC just between sender and receiver

                % loop over sender/receiver combinations to compute time-reversed GC
                iind = 0;
                for ii = 1:nchan
                    for jj = (ii+1):nchan
                        iind = iind + 1;

                        if verbose; disp(['bootstrap run ' num2str(iboot) '/' num2str(nboot) ', testing connection ' num2str(iind) '/' num2str(ninds) ': [' num2str(ii) '] -> [' num2str(jj) ']']); end

                        subset = [ii jj];

                        % autocovariance to full forward VAR model
                        [A, SIG] = autocov_to_var3(G(subset, subset, :));

                        % forward VAR model to state space VARMA models
                        [eA, eC, eK, eV, eVy] = varma2iss(reshape(A, 2, []), [], SIG, eye(2));

                        if ~isempty(intersect(output, {'TRGC', 'TRPDC', 'TRDTF'}))
                            % backward autocovariance to full backward VAR model
                            [AR, SIGR] = autocov_to_var3(permute(G(subset, subset, :), [2 1 3]));

                            % backward VAR to VARMA
                            [eAR, eCR, eKR, eVR, eVyR] = varma2iss(reshape(AR, 2, []), [], SIGR, eye(2));
                        end

                        if ~isempty(intersect(output, {'TRGC', 'GC'}))
                            % GC and TRGC computation
                            GC(:, ii, jj, iboot) = iss_SGC(eA, eC, eK, eV, z, 1, 2);

                            if ~isempty(intersect(output, {'TRGC'}))
                                GCR = iss_SGC(eAR, eCR, eKR, eVR, z, 1, 2);
                                TRGC(:, ii, jj, iboot) = GC(:, ii, jj, iboot) - GCR';
                            end

                            iind = iind + 1;
                            if verbose; disp(['bootstrap run ' num2str(iboot) '/' num2str(nboot) ', testing connection ' num2str(iind) '/' num2str(ninds) ': [' num2str(jj) '] -> [' num2str(ii) ']']); end

                            % GC and TRGC computation
                            GC(:, jj, ii, iboot) = iss_SGC(eA, eC, eK, eV, z, 2, 1);

                            if ~isempty(intersect(output, {'TRGC'}))
                                GCR = iss_SGC(eAR, eCR, eKR, eVR, z, 2, 1);
                                TRGC(:, jj, ii, iboot) = GC(:, jj, ii, iboot) - GCR';
                            end

                        end

                        % PDC
                        if ~isempty(intersect(output, {'TRPDC', 'PDC'}))
                            PDC_ = abs(est_PDC(A, freqs/2));
                            PDC(:, ii, jj, iboot) = PDC_(1, 2, :);
                            PDC(:, jj, ii, iboot) = PDC_(2, 1, :);
                        end

                        % TRPDC
                        if ~isempty(intersect(output, {'TRPDC'}))
                            TRPDC_ = PDC_ - abs(est_PDC(AR, freqs/2));
                            TRPDC(:, ii, jj, iboot) = TRPDC_(1, 2, :);
                            TRPDC(:, jj, ii, iboot) = TRPDC_(2, 1, :);
                        end

                        % DTF
                        if ~isempty(intersect(output, {'TRDTF', 'DTF'}))
                            DTF_ = abs(est_DTF(A, SIG, freqs/2));
                            DTF(:, ii, jj, iboot) = DTF_(1, 2, :);
                            DTF(:, jj, ii, iboot) = DTF_(2, 1, :);
                        end

                        % TRDTF
                        if ~isempty(intersect(output, {'TRDTF'}))
                            TRDTF_ = DTF_ - abs(est_DTF(AR, SIGR, freqs/2));
                            TRDTF(:, ii, jj, iboot) = TRDTF_(1, 2, :);
                            TRDTF(:, jj, ii, iboot) = TRDTF_(2, 1, :);
                        end

                    end
                end
            end
        end
    end
end

if ~isempty(intersect(output, {'CS'}))
    CS = permute(CS, [3 1 2 4]);
end

if ~isempty(intersect(output, {'cCOH', 'aCOH', 'iCOH'}))
    cCOH = permute(cCOH, [3 1 2 4]);
    aCOH = permute(aCOH, [3 1 2 4]);
    iCOH = permute(iCOH, [3 1 2 4]);
end

clear out
for iout = 1:length(output)
    eval(['conn.' output{iout} ' = ' output{iout} ';'])
end
conn.nlags = nlags;
