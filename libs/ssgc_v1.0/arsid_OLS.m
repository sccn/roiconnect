function [ARA,V,e] = arsid_OLS(y,p)

% Estimate a vector autoregressive model from an empirical observation time
% series using Ordinary Least Squares.
%
% y   - observation process time series
% p   - AR model order (number of lags)
%
% ARA - estimated AR coefficients array
% V   - estimated resdiuals covariance matrix
% e   - estimated residuals time series
%
% Time series y,e are matrices of dimensions n x (number of time steps), where n
% is the dimension of the observation process. The AR coefficients array ARA has
% dimensions n x n x p, so that ARA(:,:,k) is the k-lag AR coefficients matrix,
% while V is n x n symmetric and positive-(semi-)definite.

[n,T] = size(y);

assert(p < T,'model order too high (or not enough data)');

V = [];
e = [];

y = bsxfun(@minus,y,mean(y,2)); % subtract temporal mean

p1 = p+1;
Tp = T-p;

y0 = y(:,p1:T);             % lag zero
Yp = zeros(n,p,Tp);         % lags 1 through p
for k = 1:p
    Yp(:,k,:) = y(:,p1-k:T-k);
end
Yp = reshape(Yp,n*p,Tp);    % stack lagged data

ARA = y0/Yp;                % OLS ('rdivide' will generally using QR decomposition)

% IMPORTANT: test for failed OLS with something like:
%
% assert(all(isfinite(ARA(:))),'OLS failed');

if nargout > 1
    e = y0-ARA*Yp;          % residuals
    V = (e*e')/(Tp-1);      % residuals covariance matrix
    if nargout > 2
        e = [nan(n,p) e]; % pad with NaNs to align with y
    end
end

ARA = reshape(ARA,n,n,p);   % so ARA(:,:,k) is the k-lag AR coefficients matrix
