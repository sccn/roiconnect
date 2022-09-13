function PV = parcovar(V,i,k)

% Compute partial covariance PV = V(i,i) - V(i,k) * V(k,k)^{-1} * V(k,i); used
% in spectral GC calculations (see eqs. 11 and 12-15 in the reference article).
%
% V   - a (symmetric, positive-definite) covariance matrix
% i,k - non-overlapping multi-indices (vectors of integers)
%
% PV  - partial covariance matrix

if isscalar(k)
    W = sqrt(V(k,k))\V(k,i); % for 1-dim target no matrix division required
else
    opts.LT = true;                                 % we can use the more efficient linsolve (rather than rdivide)
    W = linsolve(chol(V(k,k),'lower'),V(k,i),opts); % here, since we know the Cholesky factor is lower triangular
end
PV = V(i,i)-W'*W;
