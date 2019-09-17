%% bootstrap_tsdata_to_trspwcgc
%
% Calculate sampling distribution for time-reversal-corrected 
% pairwise-conditional frequency-domain MVGCs
% from time series data, based on a nonparametric bootstrap
%
% <matlab:open('bootstrap_tsdata_to_pwcgc.m') code>
%
%% Syntax
%
%     FB = bootstrap_tsdata_to_pwcgc(U,p,nsamps,acmaxlags,acdectol)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     U          multi-trial time series data
%     p          model order (number of lags)
%     nsamps     number of bootstrap samples
%     acmaxlags  maximum autocovariance lags  (default as for 'var_to_autocov')
%     acdectol   autocovariance decay tolerance (default as for 'var_to_autocov')
%
% _output_
%
%     FB         bootstrap Granger causalities (empirical distribution)
%
%% Description
%
% Returns |nsamps| samples from the empirical sampling distribution of the
% pairwise-conditional time-domain MVGCs from the time series data |U|. The
% bootstrap randomly samples (with replacement) residuals of the full
% autoregression of |U| on its own lags; the subsampled residuals are then added
% back to the corresponding predictors to form surrogate time series [2]. |p| is
% the model order; for other parameters see <var_to_autocov.html
% |var_to_autocov|>.
%
% The first dimension of the returned matrix |FB| indexes samples, the second
% indexes the target (causee) variable and the third the source (causal)
% variable.
%
% For usage in construction of GC confidence intervals, see
% <mvgc_demo_bootstrap.html |mvgc_demo_bootstrap|>.
%
%% References
%
% [1] L. Barnett and A. K. Seth,
% <http://www.sciencedirect.com/science/article/pii/S0165027013003701 The MVGC
%     Multivariate Granger Causality Toolbox: A New Approach to Granger-causal
% Inference>, _J. Neurosci. Methods_ 223, 2014
% [ <matlab:open('mvgc_preprint.pdf') preprint> ].
%
% [2] D. A. Freedman, Bootstrapping regression models, _Ann. Stats._ 9(6), 1981.
%
%% See also
%
% <mvgc_demo_bootstrap.html |mvgc_demo_bootstrap|> |
% <bootstrap_tsdata_to_mvgc.html |bootstrap_tsdata_to_mvgc|> |
% <bootstrap_tsdata_to_smvgc.html |bootstrap_tsdata_to_smvgc|> |
% <bootstrap_tsdata_to_spwcgc.html |bootstrap_tsdata_to_spwcgc|> |
% <var_to_autocov.html |var_to_autocov|> |
% <autocov_to_pwcgc.html |autocov_to_pwcgc|>.
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function [FB, FRB, freqs] = bootstrap_tsdata_to_trspwcgc(U,p,nsamps,fres,acmaxlags,acdectol)

if nargin < 4, fres = 100; end % ensure default
if nargin < 5, acmaxlags = 2*fres-1; end % ensure default
if nargin < 6, acdectol  = []; end % ensure default

[n,m,N] = size(U);
assert(m > p,'too many lags');
p1 = p+1;
M = N*(m-p);
np = n*p;


% estimate VAR coefficients

U = demean(U);                 % no constant term

%forward
U0 = reshape(U(:,p1:m,:),n,M); % concatenate trials for unlagged observations
UL = zeros(n,p,M);
for k = 1:p
    UL(:,k,:) = reshape(U(:,p1-k:m-k,:),n,M); % concatenate trials for k-lagged observations
end
UL = reshape(UL,np,M);         % stack lags
A = U0/UL;                     % OLS using QR decomposition
if isbad(A), return; end       % something went badly wrong

%reverse
UR = flipdim(U, 2);

UR0 = reshape(UR(:,p1:m,:),n,M); % concatenate trials for unlagged observations
URL = zeros(n,p,M);
for k = 1:p
    URL(:,k,:) = reshape(UR(:,p1-k:m-k,:),n,M); % concatenate trials for k-lagged observations
end
URL = reshape(URL,np,M);         % stack lags
AR = UR0/URL;                     % OLS using QR decomposition
if isbad(AR), return; end       % something went badly wrong

% calculate predictions and residuals

m   = m-p;                     % we lose p observations
UP  = A*UL;                    % predictions
E   = U0-UP;                   % residuals: so U0 = UP + E
E   = reshape(E,n,m,N);        % put residuals back into per-trial form


URP  = AR*URL;                    % predictions
ER   = UR0-URP;                   % residuals: so U0 = UP + E
ER   = reshape(ER,n,m,N);        % put residuals back into per-trial form


FB = nan(nsamps,n,n,fres+1);

EB = zeros(n,m,N);
ERB = zeros(n,m,N);
for s = 1:nsamps
%     fprintf('PWCGC: bootstrap sample %d of %d',s,nsamps);
    
    % generate bootstrap time series
    
    in = zeros(N, m);
    for r = 1:N
        in(r, :) = randi(m,1,m);
        EB(:,:,r) = E(:,in(r, :),r); % per-trial subsample residuals with replacement
    end
    UB = UP + reshape(EB,n,M); % the bootstrap: add subsampled residuals to predictions

    inR = m-in+1;
    for r = 1:N
      ERB(:,:,r) = ER(:,inR(r, :),r); % per-trial subsample residuals with replacement
    end
    URB = URP + reshape(ERB,n,M); % the bootstrap: add subsampled residuals to predictions
    
    
    % estimate bootstrap VAR parameters
    
    AB = UB/UL;                % OLS using QR decomposition
    if isbad(AB), fprintf(' *** VAR estimation failed\n'); continue; end % something went badly wrong
    EE = UB-AB*UL;             % bootstrap residuals (recalculate)
    AB = reshape(AB,n,n,p);    % bootstrap VAR coefficients
    SIGB = (EE*EE')/(M-1);     % bootstrap residuals covariance matrix
 

    ARB = URB/URL;                % OLS using QR decomposition
    if isbad(ARB), fprintf(' *** VAR estimation failed\n'); continue; end % something went badly wrong
    EER = URB-ARB*URL;             % bootstrap residuals (recalculate)
    ARB = reshape(ARB,n,n,p);    % bootstrap VAR coefficients
    SIGRB = (EER*EER')/(M-1);     % bootstrap residuals covariance matrix
    
    
    % calculate bootstrap PWCGC
    
    [G,res] = var_to_autocov(AB,SIGB,acmaxlags,acdectol);
    if res.error, fprintf(' *** bad VAR: %s\n',res.errmsg); continue; end
    
    FB(s,:,:,:) = autocov_to_spwcgc(G, fres);
    
    
    [GR,resR] = var_to_autocov(ARB,SIGRB,acmaxlags,acdectol);
    if resR.error, fprintf(' *** bad VAR: %s\n',resR.errmsg); continue; end
    
    FRB(s,:,:,:) = autocov_to_spwcgc(GR, fres);
    
%     fprintf('\n');
end

