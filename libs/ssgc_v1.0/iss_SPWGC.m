function f = iss_SPWGC(A,C,K,V,z)

% Compute pairwise-conditional frequency-domain GCs for a state space model from
% innovations form parameters (eqs. 11 and 12-15 in the reference article).
%
% A,C,K,V - innovations form state space parameters
% z       - a vector of points on the unit circle in the complex plane
%
% f       - pairwise-conditional spectral GCs ("frequency-domain causal graph")
%
% Calculates the spectral GCs f(z,i1,i2) from y_i2 -> y_i1 conditional on y_i3,
% where i3 is the multi-index of the remaining variables in the model, for all
% pairs of (scalar) indices i1, i2. Note that for fixed z, f will not in general
% be symmetric, and will have NaNs on the diagonal.

[m,m1]  = size(A); assert(m1 == m);
[n,m1]  = size(C); assert(m1 == m);
[m1,n1] = size(K); assert(n1 == n && m1 == m);
[n1,n2] = size(V); assert(n1 == n && n2 == n);
z = z(:);

h = length(z);
f = nan(n,n,h);

H      = iss_MA(A,C,K,z);
KVSQRT = K*chol(V,'lower');
KVK    = KVSQRT*KVSQRT';

% for efficiency (in time if not memory) we pre-compute partial covariances
PV = zeros(n-1,n-1,n);
for i1 = 1:n
    oi1 = [1:i1-1 i1+1:n]; % omit i1
    PV(:,:,i1) = parcovar(V,oi1,i1);
end

for i2 = 1:n
    oi2 = [1:i2-1 i2+1:n]; % omit i2

    CR = C(oi2,:);
    [KR,VR,rep] = ss2iss(A,CR,KVK,V(oi2,oi2),K*V(:,oi2));  % reduced model Kalman gain and innovations covariance
    if rep < 0
        if     rep == -1, warning('DARE: eigenvalues on/near unit circle for source node %d',i2);
        elseif rep == -2, warning('DARE: couldn''t find stablising solution for source node %d',i2);
        end
        continue
    end
    if rep > sqrt(eps), warning('DARE: there were accuracy issues (relative residual = %e) for source node %d',rep,i2);
        continue
    end
    BR = iss_AR(A,CR,KR,z);

    for i1R = 1:n-1;
        i1  = oi2(i1R);
        oi1 = [1:i1-1 i1+1:n];  % omit i1

        SR  = VR(i1R,i1R); % reduced model spectrum is flat!
        LSR = log(SR);

        for k = 1:h
            HR = BR(i1R,:,k)*H(oi2,oi1,k); % eqs. 13, 15
            f(i1,i2,k) = LSR - log(SR-real(HR*PV(:,:,i1)*HR')); % eq. 14
        end
    end
end
