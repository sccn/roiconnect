function [paic,pbic,aic,bic] = ar_IC(y,pmax,disp)

% Compute Akaike (AIC) and Schwarz' Bayesian (BIC) information criteria for
% autoregressive model order estimation.
%
% y          - observation process time series
% pmax       - maximum AR model order (number of lags)
%
% paic, pbic - Akaike and Bayesian optimal AR model orders based on data in y
% aic,  bic  - vectors of the actual information criteria (for plotting)
%
% The time series y is a matrix of dimensions n x (number of time steps), where
% n is the dimension of the observation process. The algorithm estimates the
% likelihood by OLS very efficiently; a single QR decomposition is performed,
% and used to calculate regressions at all model orders. Any failed regressions
% will result in NaNs in the information criteria at that model order.

if nargin < 3 || isempty(disp), disp = true; end

[n,T] = size(y);

y = bsxfun(@minus,y,mean(y,2)); % subtract temporal mean

T = T-pmax; % we'll lose pmax observations

% store lags

y0 = y(:,pmax+1:pmax+T);   % lag zero
Yp = zeros(n,pmax,T);      % lags 1 through pmax
for p = 1:pmax
    Yp(:,p,:) = y(:,pmax+1-p:pmax+T-p);
end
Yp = reshape(Yp,n*pmax,T); % stack lagged data

% perform QR decomposition for all lags in one shot

[Q,R] = qr(Yp',0); % "economy-size" decomposition is all we need
y0Q = y0*Q;

aic = nan(pmax,1);
bic = nan(pmax,1);

% loop through model orders

for p = 1:pmax
    np = n*p;
    r = min(np,T);
    ARA = y0Q(:,1:r)/R(1:r,1:np)';      % OLS for regression against first p lags only
    if ~all(isfinite(ARA)), continue; end % show-stopper - skip and carry on
    e = y0-ARA*Yp(1:np,:);              % residuals
    L = log(det((e*e')/(T-1)));         % likelihood
    nfp = n*np;                         % number of free parameters
    if T > nfp+1
        aic(p) = L + 2*nfp/(T-nfp-1);   % with Hurvich and Tsai small-sample correction
    end
    bic(p) = L + nfp*log(T)/T;
end

% optimal model orders (note: NaNs are ignored)

morder = 1:pmax;
[~,idx] = min(aic); paic = morder(idx);
[~,idx] = min(bic); pbic = morder(idx);

if disp
    plot(morder,[aic bic]);
    xlim([1 pmax]);
    title('AR model order information criteria');
    legend(sprintf('AIC (opt = %d)',paic),sprintf('BIC (opt = %d)',pbic),'location','northwest');
end
