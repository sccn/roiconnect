function F = iss_GC(A,C,K,V,i1,i2)

% Compute time-domain (conditional) GC for a state space model from innovations
% form parameters (eq. 9 in the reference article).
%
% A,C,K,V - innovations form state space parameters
% i1      - target variable multi-index (vector of indices)
% i2      - source variable multi-index (vector of indices)
%
% F       - GC from i2 -> i1
%
% Calculates GC from y_i2 -> y_i1 conditional on y_i3, where i3 is the
% multi-index of the remaining variables in the model. The multi-indices i1, i2
% must not overlap.

[m,m1]  = size(A); assert(m1 == m);
[n,m1]  = size(C); assert(m1 == m);
[m1,n1] = size(K); assert(n1 == n && m1 == m);
[n1,n2] = size(V); assert(n1 == n && n2 == n);

i1 = i1(:)'; assert(all(i1 >=1 & i1 <= n));
i2 = i2(:)'; assert(all(i2 >=1 & i2 <= n));
assert(isempty(intersect(i1,i2)));

i3  = 1:n; i3([i1 i2]) = []; % indices of remaining (conditioning) variables
i13 = [i1 i3];
i1R = 1:length(i1);

F = NaN;

KVSQRT = K*chol(V,'lower');
[~,VR,rep] = ss2iss(A,C(i13,:),KVSQRT*KVSQRT',V(i13,i13),K*V(:,i13)); % reduced model innovations covariance
if rep < 0
    if     rep == -1, warning('DARE: eigenvalues on/near unit circle');
    elseif rep == -2, warning('DARE: couldn''t find stablising solution');
    end
    return
end
if rep > sqrt(eps), warning('DARE: there were accuracy issues (relative residual = %e)',rep);
    return
end

F = log(det(VR(i1R,i1R))) - log(det(V(i1,i1))); % eq. 9
