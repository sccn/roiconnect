function [y,e] = ar_gen(ARA,V,T)

% Generate a realization of a vector autoregressive process with Gaussian
% residuals.
%
% ARA - AR coefficients array
% V   - residuals covariance matrix
% T   - sequence length (time steps)
%
% y   - observation time series
% e   - residuals time series
%
% Time series y,e are matrices of dimensions n x (number of time steps), where n
% is the dimension of the observation process. The AR coefficients array ARA has
% dimensions n x n x p, so that ARA(:,:,k) is the k-lag AR coefficients matrix,
% while V is n x n.
%
% Note that we need rho < 1 for a stable process, where rho is the spectral norm
% of ARA (see routine specnorm.m). V must be symmetric positive-definite.
%
% The returned time series should be truncated at the beginning if an
% (approximately) stationary process is required. Since transients for a stable
% AR(p) process decay exponentially at a rate rho, a reasonable number of time
% steps to truncate is of the order log(eps)/log(rho).

[n,n1,p] = size(ARA); assert(n1 == n);

% Gaussian white noise innovations with covariance matrix V

e = chol(V,'lower')*randn(n,T);

% Construct observations - autoregress

y = e;
for t = p+1:T
    for k = 1:p
        y(:,t) = y(:,t) + ARA(:,:,k)*y(:,t-k);
    end
end
