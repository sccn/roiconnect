function H = iss_MA(A,C,K,z)

% Compute moving-average representation (transfer function) for a state space
% model in innovations form (eq. 4 in the reference article).
%
% A,C,K - innovations form state space parameters
% z     - a vector of points on the unit circle in the complex plane
%
% H     - transfer function

[m,m1]  = size(A); assert(m1 == m);
[n,m1]  = size(C); assert(m1 == m);
[m1,n1] = size(K); assert(n1 == n && m1 == m);
z = z(:);

h = length(z);
In = eye(n);
Im = eye(m);
H = zeros(n,n,h);
for k = 1:h
    H(:,:,k) = In + C*((z(k)*Im-A)\K); % eq. 4
end
