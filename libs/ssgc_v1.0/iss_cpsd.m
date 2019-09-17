function S = iss_cpsd(H,V)

% Calculate cross-power spectral density from transfer function and innovations
% covariance matrix (eq. 6 in the reference article).
%
% H - transfer function
% V - innovations covariance matrix
%
% S - cross-power spectral density array (last dimension is frequency)
%
% You may calculate the cpsd from innovations form parameters A,C,K,V by calling
%
% S = iss_cpsd(iss_MA(A,C,K,z),V);

[n,n1,h] = size(H); assert(n1 == n);
[n1,n2]  = size(V); assert(n1 == n && n2 == n);

S = zeros(n,n,h);
SQRTV = chol(V,'lower');
for k = 1:h
    Hk = H(:,:,k)*SQRTV;
    S(:,:,k) = Hk*Hk';
end
