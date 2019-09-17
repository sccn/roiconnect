function [m,A,C,K,V,z,e] = s4sid_CCA(y,pf,mo,disp)

% Estimate an (innovations form) state space model from an empirical observation
% time series using Larimore's Canonical Correlations Analysis (CCA) state
% space-subspace algorithm (W. E. Larimore, in Proc. Amer. Control Conference,
% vol. 2, IEEE, 1983).
%
% y       - observation process time series
% pf      - past/future horizons for canonical correlations
% mo      - model order (see note below)
% disp    - flag: set to true (default) to display singular values and SVC
%
% m       - model order deployed
% A,C,K,V - estimated innovations form state space parameters
% z       - estimated state process time series
% e       - estimated innovations process time series
%
% Time series y,z,e are matrices of dimensions n x (number of time steps), where
% n is the dimension of the observation process.
%
% The past/future horizons pf may be supplied as a 2-vector [p,f] or a scaler p
% = f = pf. Bauer recommends setting p = f = 2*pAIC, where pAIC is the optimal
% AR model order for the observation process according to Aikaike's Information
% Criterion, although we have obtained good results using p = f = pBIC, where
% pBIC is the optimal AR model order according to Schwarz' Bayesian Information
% Criterion.
%
% If the model order mo is empty (or unspecified), Bauer's SVC model order
% selection criterion (D. Bauer, Automatica 37, 1561, 2001) is used to estimate
% an optimal model order m. If mo is set to zero, the CCA singular values are
% displayed and the user is prompted to enter an (integer) model order.
% Otherwise, mo may be set to any positive integer model order <= n*min(p,f),
% where n is the number of observation variables.

[n,T] = size(y);

if isscalar(pf);
    p = pf;    f = pf;
elseif isvector(pf) && length(pf) == 2
    p = pf(1); f = pf(2);
else
    error('past/future horizon must be a 2-vector or a scalar');
end
assert(p+f < T,'past/future horizon too large (or not enough data)');
momax = n*min(p,f);

if nargin < 3,                  mo   = [];   end % default: use Bauer's SVC model order estimate
if nargin < 4 || isempty(disp), disp = true; end % default: display singular values 

assert(isempty(mo) || isequal(mo,0) || (isscalar(mo) && mo == floor(mo) && mo > 0 && mo <= momax),'model order must be empty (use SVC), 0 (prompt), or a positive integer <= n*min(p,f) = %d',momax);

prompt = isequal(mo,0); % display singular values and prompt user for model order
useSVC = isempty(mo);   % use Bauer's SVC model order estimate
if prompt, disp = true; end

A = [];
C = [];
K = [];
V = [];
z = [];
e = [];

y = bsxfun(@minus,y,mean(y,2)); % subtract temporal mean

Tp  = T-p+1;
Tf  = T-f;
Tpf = Tp-f;

Yf = zeros(n,f,Tpf);
for k = 1:f
    Yf(:,k,:) = y(:,p+k:Tf+k);
end
Yf = reshape(Yf,n*f,Tpf);

YP = zeros(n,p,Tp);
for k = 0:p-1
    YP(:,k+1,:) = y(:,p-k:T-k);
end
YP = reshape(YP,n*p,Tp);

Yp = YP(:,1:Tpf);

[Wf,cholp] = chol((Yf*Yf')/Tpf,'lower');
assert(cholp == 0, 'forward weight matrix not positive-definite');

[Wp,cholp] = chol((Yp*Yp')/Tpf,'lower');
assert(cholp == 0, 'backward weight matrix not positive-definite');

BETA = Yf/Yp; % 'OH' estimate: regress future on past
assert(all(isfinite(BETA(:))),'subspace regression failed');

[~,S,U] = svd(Wf\BETA*Wp); % SVD of CCA-weighted OH estimate

sval = diag(S); % the singular values

if disp || useSVC
    nparms = 2*n*(1:momax)';                           % number of free parameters (Hannan & Deistler, see also Bauer 2001)
    SVC = -log(1-[sval(2:end);0]) + nparms*(log(T)/T); % Bauer's Singular Value Criterion
    [~,moSVC] = min(SVC);
end

if disp
    display_mo(SVC,moSVC,sval)
end

if prompt
    m = prompt_mo(moSVC,momax);
elseif useSVC
    m = moSVC;
else
    m = mo;
end

if nargout < 2, return; end % just want model order

z = (diag(sqrt(sval(1:m)))*U(:,1:m)'/Wp)*YP; % Kalman states estimate

% Calculate model parameters by regression
    
C = y(:,p+1:T)/z(:,1:Tp-1);     % observation matrix
assert(all(isfinite(C(:))), 'C parameter estimation failed');

e = y(:,p+1:T) - C*z(:,1:Tp-1); % innovations
V = (e*e')/(Tp-2);              % innovations covariance matrix (unbiased estimate)

AK = z(:,2:Tp)/[z(:,1:Tp-1);e];
assert(all(isfinite(AK(:))), 'A,K parameter estimation failed');

A = AK(:,1:m);                  % state transition matrix
K = AK(:,m+1:m+n);              % Kalman gain matrix

end

function display_mo(SVC,moSVC,sval)

    monum = length(sval);
    mos = (1:monum)';
    
    subplot(2,1,1); % singular values
    bar(mos,sval,'b');
    xlim([0 monum+1]);
    xlabel('model order');
    ylabel('singular value');
    title('singular values');
    
    subplot(2,1,2); % SVC
    icmin = min(SVC);
    icmax = max(SVC);
    SVC = (SVC-icmin)/(icmax-icmin); % scale between 0 and 1
    plot(mos,SVC);
    xlim([0 monum+1]);
    xlabel('model order');
    ylabel('information criterion');
    title(sprintf('SVC (optimum model order = %d)',moSVC));

end

function m = prompt_mo(moSVC,momax)

    while true
        mo_user = input(sprintf('RETURN for SVC model order = %d, or enter a numerical model order <= %d: ',moSVC,momax),'s');
        if isempty(mo_user), m = moSVC; return; end   
        m = str2double(mo_user);
        if ~isnan(m) && isscalar(m) && floor(m) == m && m > 0 && m <= momax, return; end
        fprintf(2,'ERROR: ');
    end

end
