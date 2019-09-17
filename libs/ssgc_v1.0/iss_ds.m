function [Ak,Kk,Vk] = iss_ds(A,C,K,V,k)

% Calculate innovations form state space parameters for a downsampled state
% space model
%
% A,C,K,V  - innovations form state space parameters
% k        - downsample factor (positive integer)
%
% Ak,Kk,Vk - innovations form state space parameters for downsampled model
%
% See V. Solo, arXiv:1501.04663v1, Sec. 4.

[m,m1]  = size(A); assert(m1 == m);
[n,m1]  = size(C); assert(m1 == m);
[m1,n1] = size(K); assert(n1 == n && m1 == m);
[n1,n2] = size(V); assert(n1 == n && n2 == n);
assert(isscalar(k) && isnumeric(k) && k == floor(k) && k > 0);

if k == 1 % nothing to do
    Ak = A;
    Kk = K;
    Vk = V;
    return
end

Ak = eye(m);   % Ak = A^0
AkKSQRTV = K*chol(V,'lower');
Qk = AkKSQRTV*AkKSQRTV';
for p = 1:k-1
    Ak = Ak*A; % Ak = A^p
    AkKSQRTV = Ak*AkKSQRTV;
    Qk = Qk + AkKSQRTV*AkKSQRTV';
end
Sk = Ak*K*V;   % Ak = A^(k-1)
Ak = Ak*A;     % Ak = A^k

% We now have general form state space parameters for the downsampled model
% (note: Ck = C and Rk = V for any k). We solve the DARE to get remaining
% innovations form parameters.

[Kk,Vk] = ss2iss(Ak,C,Qk,V,Sk);
