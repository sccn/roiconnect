function conn = data2sctrgcmim(data, fres, nlags, cond, nboot, maxfreq, inds, output, verbose)
% Epoched time series data to spectral conditional time-reversed Granger causality
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
% influence of weak data asymmetries on Granger-causal analyses. 
% In International Conference on Latent Variable Analysis and Signal 
% Separation (pp. 25-33). Springer, Berlin, Heidelberg.
%
% Winkler, I., Panknin, D., Bartz, D., M?ller, K. R., & Haufe, S. (2016). 
% Validity of time reversal for testing Granger causality. IEEE Transactions 
% on Signal Processing, 64(11), 2746-2760.
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


[nchan, ndat, nepo] = size(data);

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

freqs = freqs(1:maxfreq);
z = exp(-i*pi*freqs);

if nargin < 7 || isempty(inds)
  inds = {}; ninds = 0;
  for ii = 1:nchan
    for ij = (ii+1):nchan
      inds{ninds+1} = {ii, ij};    
%       inds{ninds+2} = {ij, ii};  
      ninds = ninds + 1;
    end
  end
else
  ninds = length(inds);  
end

if nargin < 8 || isempty(output)
  output = {'TRGC'};
end

if nargin < 9 || isempty(verbose)
  verbose = 0;
end

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
  

clear TRGC GC MIM MIC CS COH wPLI

if abs(nboot) < 1 % no bootstrap

  % data to autocovariance
  %     G = tsdata_to_autocov(data, nlags);
  [CS, ~, wPLI, ~] = data2cs_event(data(:, :)', ndat, floor(ndat/2), ndat, [], CSpara);

  maxfreq = size(CS,3);
      
  if ~isempty(intersect(output, {'MIM', 'MIC', 'COH'}))
    clear COH
    for ifreq = 1:maxfreq
        clear pow
        pow = real(diag(CS(:,:,ifreq)));
        COH(:,:,ifreq) = CS(:,:,ifreq)./ sqrt(pow*pow');
    end
  end
  
  if ~isempty(intersect(output, {'GC', 'TRGC'}))
    G = cpsd_to_autocov(CS, nlags);
  end
  
  if isempty(intersect(output, {'CS'}))
    clear CS
  end

  if cond 
  % (time-reversed) GC conditioned on all other variables

  
    if ~isempty(intersect(output, {'GC', 'TRGC'}))
      % autocovariance to full forward VAR model
      [A, SIG] = autocov_to_var3(G);

      % forward VAR model to state space VARMA models
      [eA2, eC2, eK2, eV2, eVy2] = varma2iss(reshape(A, nchan, []), [], SIG, eye(nchan)); 

      if ~isempty(intersect(output, {'TRGC'}))
        % backward autocovariance to full backward VAR model
        [AR, SIGR] = autocov_to_var3(permute(G, [2 1 3]));

        % backward VAR to VARMA
        [eA2R, eC2R, eK2R, eV2R, eVy2R] = varma2iss(reshape(AR, nchan, []), [], SIGR, eye(nchan));
      end
    end

    %% loop over sender/receiver combinations to compute (time-reversed) GC 
    for iind = 1:ninds   
      if ~isequal(inds{iind}{1}, inds{iind}{2})
        if verbose
          disp(['testing connection ' num2str(iind) '/' num2str(ninds) ': [' num2str(inds{iind}{1}) '] <-> [' num2str(inds{iind}{2}) '], conditional'])
        end
        
        if ~isempty(intersect(output, {'MIM', 'MIC'}))
        %MIC and MIM
          [MIC(:, iind) , MIM(:, iind)] =  roi_mim2(COH, inds{iind}{1}, inds{iind}{2});    
        end
        
        if ~isempty(intersect(output, {'GC', 'TRGC'}))
          GC(:, iind, 1) = iss_SGC(eA2, eC2, eK2, eV2, z, inds{iind}{2}, inds{iind}{1});
          GC(:, iind, 2) = iss_SGC(eA2, eC2, eK2, eV2, z, inds{iind}{1}, inds{iind}{2});
          
          if ~isempty(intersect(output, {'TRGC'}))
            GCR = iss_SGC(eA2R, eC2R, eK2R, eV2R, z, inds{iind}{2}, inds{iind}{1});
            TRGC(:, iind, 1) = GC(:, iind, 1) - GCR';
            GCR = iss_SGC(eA2R, eC2R, eK2R, eV2R, z, inds{iind}{1}, inds{iind}{2});
            TRGC(:, iind, 2) = GC(:, iind, 2) - GCR';
          end
        end
      else       
        if ~isempty(intersect(output, {'GC', 'TRGC'}))
          GC(:, iind, 1:2) = 0;
          TRGC(:, iind, 1:2) = 0;
        end
      end
    end      
  else
  % (time-reversed) GC just between sender and receiver sets

    % loop over sender/receiver combinations to compute time-reversed GC 
    fprintf('Progress of %d:', ninds);
    for iind = 1:ninds  
        if mod(iind,100) == 0
            fprintf('%d', iind);
        elseif mod(iind, 10) == 0
            fprintf('.');
        end
      if ~isequal(inds{iind}{1}, inds{iind}{2})
        if verbose
          disp(['testing connection ' num2str(iind) '/' num2str(ninds) ': [' num2str(inds{iind}{1}) '] <-> [' num2str(inds{iind}{2}) ']'])
        end
        %ind configuration 
        subset = [inds{iind}{1} inds{iind}{2}];
        nsubsetvars = length(subset);
        subinds = {1:length(inds{iind}{1}), length(inds{iind}{1}) + (1:length(inds{iind}{2}))};
     
        if ~isempty(intersect(output, {'MIM', 'MIC'}))
          %MIC and MIM
          [MIC(:, iind) , MIM(:, iind)] =  roi_mim2(COH(subset, subset, :), subinds{1}, subinds{2});
        end
        
        if ~isempty(intersect(output, {'GC', 'TRGC'}))
          % autocovariance to full forward VAR model
          [A, SIG] = autocov_to_var3(G(subset, subset, :));

          % forward VAR model to state space VARMA models
          [eA2, eC2, eK2, eV2, eVy2] = varma2iss(reshape(A, nsubsetvars, []), [], SIG, eye(nsubsetvars)); 

          % GC and TRGC computation
          GC(:, iind, 1) = iss_SGC(eA2, eC2, eK2, eV2, z, subinds{2}, subinds{1});
          GC(:, iind, 2) = iss_SGC(eA2, eC2, eK2, eV2, z, subinds{1}, subinds{2});
          
          if ~isempty(intersect(output, {'TRGC'}))
            % backward autocovariance to full backward VAR model
            [AR, SIGR] = autocov_to_var3(permute(G(subset, subset, :), [2 1 3]));

            % backward VAR to VARMA
            [eA2R, eC2R, eK2R, eV2R, eVy2R] = varma2iss(reshape(AR, nsubsetvars, []), [], SIGR, eye(nsubsetvars));
            GCR = iss_SGC(eA2R, eC2R, eK2R, eV2R, z, subinds{2}, subinds{1});
            TRGC(:, iind, 1) = GC(:, iind, 1) - GCR';
            GCR = iss_SGC(eA2R, eC2R, eK2R, eV2R, z, subinds{1}, subinds{2});
            TRGC(:, iind, 2) = GC(:, iind, 2) - GCR';
          end
        end
%         
%         F = autocov_to_pwcgc(G(subset, subset, :));
%         FR = autocov_to_pwcgc(permute(G(subset, subset, :), [2 1 3]));
      else
        if ~isempty(intersect(output, {'GC', 'TRGC'}))
          GC(:, iind, 1:2) = 0;
          TRGC(:, iind, 1:2) = 0;
        end
      end
    end
    fprintf('\n');
  end
else % bootstrap

  % loop over bootstrap samples
  for iboot = 1:nboot
    bootinds = randi(nepo, nepo, 1);

    % data to autocovariance
%       G = tsdata_to_autocov(data(:, :, bootinds), nlags);
%     CS = tsdata_to_cpsd(data(:, :, bootinds), fres, 'MT', ndat); 
    data_ = data(:, :, bootinds);
    CS_ = data2cs_event(data_(:, :)', ndat, ndat, ndat, fres+1, CSpara);
    
    if ~isempty(intersect(output, {'CS'})) 
      CS(:, :, :, iboot) = CS_(:, :, 1:maxfreq); 
    end
    
    if ~isempty(intersect(output, {'MIM', 'MIC', 'COH'}))
      clear COH_
      for ifreq = 1:maxfreq
          clear pow
          pow = real(diag(CS_(:,:,ifreq)));
          COH_(:,:,ifreq) = CS_(:,:,ifreq)./ sqrt(pow*pow');
      end
      if ~isempty(intersect(output, {'COH'})) 
        COH(:, :, :, iboot) = COH_;
      end
    end

    if ~isempty(intersect(output, {'GC', 'TRGC'}))
      G = cpsd_to_autocov(CS_, nlags);
    end
    
    clear CS_

    if cond 
    % (time-reversed) GC conditioned on all other variables

      if ~isempty(intersect(output, {'GC', 'TRGC'}))
        % autocovariance to full forward VAR model
        [A, SIG] = autocov_to_var3(G);

        % forward VAR model to state space VARMA models
        [eA, eC, eK, eV, eVy] = varma2iss(reshape(A, nchan, []), [], SIG, eye(nchan)); 

        if ~isempty(intersect(output, {'TRGC'}))
          % backward autocovariance to full backward VAR model
          [AR, SIGR] = autocov_to_var3(permute(G, [2 1 3]));

          % backward VAR to VARMA
          [eAR, eCR, eKR, eVR, eVyR] = varma2iss(reshape(AR, nchan, []), [], SIGR, eye(nchan));
        end
      end

      %% loop over sender/receiver combinations to compute (time-reversed) GC 
      for iind = 1:ninds
        if ~isequal(inds{iind}{1}, inds{iind}{2})
%           i3 = setdiff(setdiff(1:nchan, inds{iind}{1}), inds{iind}{2});
          if verbose
            disp(['bootstrap run ' num2str(iboot) '/' num2str(nboot) ', testing connection ' num2str(iind) '/' num2str(ninds) ': [' num2str(inds{iind}{1}) '] <-> [' num2str(inds{iind}{2}) '], conditional'])
          end
          
          if ~isempty(intersect(output, {'MIM', 'MIC'}))
            %MIC and MIM
            [MIC(:, iind, iboot) , MIM(:, iind, iboot)] =  roi_mim2(COH_, inds{iind}{1}, inds{iind}{2});  
          end
          
          if ~isempty(intersect(output, {'GC', 'TRGC'}))
            GC(:, iind, 1, iboot) = iss_SGC(eA, eC, eK, eV, z, inds{iind}{2}, inds{iind}{1});
            GC(:, iind, 2, iboot) = iss_SGC(eA, eC, eK, eV, z, inds{iind}{1}, inds{iind}{2});
            
            if ~isempty(intersect(output, {'TRGC'}))
              GCR = iss_SGC(eAR, eCR, eKR, eVR, z, inds{iind}{2}, inds{iind}{1});
              TRGC(:, iind, 1, iboot) = GC(:, iind, 1, iboot) - GCR';
              GCR = iss_SGC(eAR, eCR, eKR, eVR, z, inds{iind}{1}, inds{iind}{2});
              TRGC(:, iind, 2, iboot) = GC(:, iind, 2, iboot) - GCR';
            end
          end
        else
          if ~isempty(intersect(output, {'GC', 'TRGC'}))
            GC(:, iind, 1:2, iboot) = 0;
            TRGC(:, iind, 1:2, iboot) = 0;
          end
        end
      end
    else
    % (time-reversed) GC just between sender and receiver sets

      % loop over sender/receiver combinations to compute time-reversed GC 
      for iind = 1:ninds 
        if ~isequal(inds{iind}{1}, inds{iind}{2})
          if verbose
            disp(['bootstrap run ' num2str(iboot) '/' num2str(nboot) ', testing connection ' num2str(iind) '/' num2str(ninds) ': [' num2str(inds{iind}{1}) '] <-> [' num2str(inds{iind}{2}) ']'])
          end
          
          subset = [inds{iind}{1} inds{iind}{2}];
          nsubsetvars = length(subset);
          subinds = {1:length(inds{iind}{1}), length(inds{iind}{1}) + (1:length(inds{iind}{2}))};   
          
          if ~isempty(intersect(output, {'MIM', 'MIC'}))
            %MIC and MIM
            [MIC(:, iind, iboot) , MIM(:, iind, iboot)] =  roi_mim2(COH_(subset, subset, :), subinds{1}, subinds{2});
          end

          if ~isempty(intersect(output, {'GC', 'TRGC'}))
            % autocovariance to full forward VAR model
            [A, SIG] = autocov_to_var3(G(subset, subset, :));

            % forward VAR model to state space VARMA models
            [eA, eC, eK, eV, eVy] = varma2iss(reshape(A, nsubsetvars, []), [], SIG, eye(nsubsetvars)); 

            % GC and TRGC computation
            GC(:, iind, 1, iboot) = iss_SGC(eA, eC, eK, eV, z, subinds{2}, subinds{1});
            GC(:, iind, 2, iboot) = iss_SGC(eA, eC, eK, eV, z, subinds{1}, subinds{2});
            
            if ~isempty(intersect(output, {'TRGC'}))
              % backward autocovariance to full backward VAR model
              [AR, SIGR] = autocov_to_var3(permute(G(subset, subset, :), [2 1 3]));

              % backward VAR to VARMA
              [eAR, eCR, eKR, eVR, eVyR] = varma2iss(reshape(AR, nsubsetvars, []), [], SIGR, eye(nsubsetvars));
              GCR = iss_SGC(eAR, eCR, eKR, eVR, z, subinds{2}, subinds{1});
              TRGC(:, iind, 1, iboot) = GC(:, iind, 1, iboot) - GCR';
              GCR = iss_SGC(eAR, eCR, eKR, eVR, z, subinds{1}, subinds{2});
              TRGC(:, iind, 2, iboot) = GC(:, iind, 2, iboot) - GCR';
            end
          end
        else
          if ~isempty(intersect(output, {'GC', 'TRGC'}))
            GC(:, iind, 1:2, iboot) = 0;
            TRGC(:, iind, 1:2, iboot) = 0;
          end
        end
      end        
    end       
  end    
end

if ~isempty(intersect(output, {'CS'})) 
  CS = permute(CS, [3 1 2 4]);
end

if ~isempty(intersect(output, {'COH'})) 
  COH = permute(COH, [3 1 2 4]);
end

clear out
for iout = 1:length(output)
  eval(['conn.' output{iout} ' = ' output{iout} ';'])
end
conn.nlags = nlags;
conn.inds = inds;
