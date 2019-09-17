function [y,z,e] = iss_gen(A,C,K,V,T)

% Generate a realization of a state space process with Gaussian innovations from
% innovations form state space parameters.
%
% A,C,K,V - innovations form state space parameters
% T       - sequence length (time steps)
%
% y       - observation time series
% z       - state time series
% e       - innovations time series
%
% Time series y,z,e are matrices of dimensions (number of variables) x (number
% of time steps).
%
% Note that we need rhoA < 1 for a stable process and rhoB < 1 for a minimum
% phase process, where rhoA is the spectral norm of A and rhoB the spectral norm
% of B = A-K*C (see routine specnorm.m). V must be symmetric positive-definite.
%
% The returned time series should be truncated at the beginning if an
% (approximately) stationary process is required. Since transients for the
% stable AR(1) state space process z decay exponentially at a rate rhoA, a
% reasonable number of time steps to truncate is of the order
% log(eps)/log(rhoA).

[m,m1]  = size(A); assert(m1 == m);
[n,m1]  = size(C); assert(m1 == m);
[m1,n1] = size(K); assert(n1 == n && m1 == m);
[n1,n2] = size(V); assert(n1 == n && n2 == n);

% Gaussian white noise innovations with covariance matrix V

e = chol(V,'lower')*randn(n,T);

% z is AR(1) with residuals K*e

z =  K*e;
for t = 2:T
    z(:,t) = A*z(:,t-1) + z(:,t);
end

% Construct observations

y = C*z + e;


