function B = iss_AR(A,C,K,z)

% Compute autoregressive representation (inverse transfer function) for a state
% space model in innovations form (eq. 5 in the reference article).
%
% A,C,K - innovations form state space parameters
% z     - a vector of points on the unit circle in the complex plane
%
% B     - inverse transfer function

[m,m1]  = size(A); assert(m1 == m);
[n,m1]  = size(C); assert(m1 == m);
[m1,n1] = size(K); assert(n1 == n && m1 == m);
z = z(:);

h = length(z);
In = eye(n);
Im = eye(m);
BB = A-K*C;
B = zeros(n,n,h);
for k = 1:h
    B(:,:,k) = In - C*((z(k)*Im-BB)\K); % eq. 5
end
