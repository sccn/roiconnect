function rho = specnorm(ARA)

% Calculate spectral norm (maximum absolute eigenvalue of the associated AR(1)
% transition matrix) for a vector autoregressive model.
%
% ARA - AR coefficients array
%
% rho - spectral norm of ARA
%
% ARA is n x n x p where where n is the dimension of the model and p the model
% order, so that ARA(:,:,k) is the k-lag AR coefficients matrix. The model is
% stable iff rho < 1.
%
% Note that for a state space model in innovations form with parameters (A,C,K),
% the model is stable iff specnorm(A) < 1, and minimum-phase iff specnorm(B) <
% 1, where B = A - K*C.

[n,n1,p] = size(ARA); assert(n1 == n); % n is dimensionality, p is model order

if p > 1 % set up associated AR(1) matrix
    ARA = [reshape(ARA,n,p*n); eye((p-1)*n) zeros((p-1)*n,n)];
end
rho = max(abs(eig(ARA,'nobalance')));
