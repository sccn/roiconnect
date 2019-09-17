function [A,C,K,rho] = ar2iss(ARA)
    
% Return innovations form state space model parameters corresponding to a vector
% autoregressive model.
%
% ARA   - AR coefficients array
%
% A,C,K - innovations form state space parameters
%
% rho   - AR spectral norm
%
% ARA is an n x n x p matrix, where n is the dimension of the observation
% variable and p the AR model order, so that ARA(:,:,k) is the AR coefficients
% matrix at lag k. Note that rho >= 1 indicates an unstable AR process: rho > 1
% is explosive, rho close to 1 may be unit-root.

[n,n1,p] = size(ARA); % p is VAR model order
assert(n1 == n);
pn1 = (p-1)*n;

C = reshape(ARA,n,p*n);
A = [C; eye(pn1) zeros(pn1,n)];
K = [eye(n); zeros(pn1,n)];

if nargout > 3
    rho = max(abs(eig(A,'nobalance')));
end
