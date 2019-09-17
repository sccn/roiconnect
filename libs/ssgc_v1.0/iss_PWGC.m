function F = iss_PWGC(A,C,K,V)

% Compute pairwise-conditional time-domain GCs for a state space model from
% innovations form parameters (eq. 9 in the reference article).
%
% A,C,K,V - innovations form state space parameters
%
% F       - pairwise-conditional GCs ("causal graph")
%
% Calculates the GCs F(i1,i2) from y_i2 -> y_i1 conditional on y_i3, where i3 is
% the multi-index of the remaining variables in the model, for all pairs of
% (scalar) indices i1, i2. Note that F will not in general be symmetric, and
% will have NaNs on the diagonal.

[m,m1]  = size(A); assert(m1 == m);
[n,m1]  = size(C); assert(m1 == m);
[m1,n1] = size(K); assert(n1 == n && m1 == m);
[n1,n2] = size(V); assert(n1 == n && n2 == n);

F = nan(n);

KVSQRT = K*chol(V,'lower');
KVK    = KVSQRT*KVSQRT';
LV     = log(diag(V));

for i2 = 1:n
    oi2 = [1:i2-1 i2+1:n]; % omit i2
    
    [~,VR,rep] = ss2iss(A,C(oi2,:),KVK,V(oi2,oi2),K*V(:,oi2)); % reduced model innovations covariance
    if rep < 0
        if     rep == -1, warning('DARE: eigenvalues on/near unit circle for source %d',i2);
        elseif rep == -2, warning('DARE: couldn''t find stablising solution for source node %d',i2);
        end
        continue
    end
    if rep > sqrt(eps), warning('DARE: there were accuracy issues (relative residual = %e) for source %d',rep,i2);
        continue
    end
    
    F(oi2,i2) = log(diag(VR))-LV(oi2);
end
