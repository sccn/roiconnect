function ARA = ar_rand(n,p,rhoa)

% Generate random stable normally distributed autoregression coefficients.
%
% n    - observation variable dimension
% p    - AR model order
% rhoa - spectral norm
%
% ARA  - autoregression coefficients array
%
% We need rhoa < 1 for a stable model. This routine creates ARA as a random
% Gaussian array of dimensions n x n x p and scales it so that its spectral norm
% is equal to rhoa.

assert(rhoa < 1);

ARA = randn(n,n,p);
dfac = rhoa/specnorm(ARA);
f = 1;
for k = 1:p
    f = dfac*f;
    ARA(:,:,k) = f*ARA(:,:,k);
end
