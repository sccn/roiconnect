function f = iss_SGC(A,C,K,V,z,i1,i2)

% Compute frequency-domain (conditional) GC for a state space model from innovations
% form parameters (eqs. 11 and 12-15 in the reference article).
%
% A,C,K,V - innovations form state space parameters
% z       - a vector of points on the unit circle in the complex plane
% i1      - target variable multi-index (vector of indices)
% i2      - source variable multi-index (vector of indices)
%
% f       - spectral GC from i2 -> i1
%
% Calculates spectral GC from y_i2 -> y_i1 conditional on y_i3, where i3 is the
% multi-index of the remaining variables in the model. The multi-indices i1, i2
% must not overlap.

[m,m1]  = size(A); assert(m1 == m);
[n,m1]  = size(C); assert(m1 == m);
[m1,n1] = size(K); assert(n1 == n && m1 == m);
[n1,n2] = size(V); assert(n1 == n && n2 == n);
z = z(:);

i1 = i1(:)'; assert(all(i1 >=1 & i1 <= n));
i2 = i2(:)'; assert(all(i2 >=1 & i2 <= n));
assert(isempty(intersect(i1,i2)));

i3  = 1:n; i3([i1 i2]) = [];
i13 = [i1 i3];
i23 = [i2 i3];
i1R = 1:length(i1);

h = length(z);
f = nan(1,h);

H      = iss_MA(A,C,K,z);
VSQRT  = chol(V,'lower');
PVSQRT = chol(parcovar(V,i23,i1),'lower');

% Note: in theory we shouldn't have to take the real part of the determinants of
% the (Hermitian, positive-definite) matrices in the calculation of the f(k),
% since they should be real. However, Matlab's det() will sometimes return a
% tiny imaginary part.

if isempty(i3) % unconditional (note: does not require reduced model)

    for k = 1:h
        HV   = H(:,:,k)*VSQRT;
        S    = HV*HV'; % CPSD (eq. 6)
        S11  = S(i1,i1);
        HV12 = H(i1,i2,k)*PVSQRT;
        f(k) = log(real(det(S11))) - log(real(det(S11-HV12*HV12'))); % eq. 11
    end

else % conditional

    KVSQRT   = K*VSQRT;
    CR       = C(i13,:);
    [KR,VR,rep]  = ss2iss(A,CR,KVSQRT*KVSQRT',V(i13,i13),K*V(:,i13)); % reduced model Kalman gain and innovations covariance
    if rep < 0
        if     rep == -1, warning('DARE: eigenvalues on/near unit circle');
        elseif rep == -2, warning('DARE: couldn''t find stablising solution');
        end
        return
    end
    if rep > sqrt(eps), warning('DARE: there were accuracy issues (relative residual = %e)',rep);
        return
    end
    BR       = iss_AR(A,CR,KR,z);
    SR       = VR(i1R,i1R); % reduced model spectrum is flat!
    LDSR     = log(det(SR));
    for k = 1:h
        HR   = BR(i1R,:,k)*H(i13,i23,k)*PVSQRT;  % eqs. 13, 15
        f(k) = LDSR - log(real(det(SR-HR*HR'))); % eq. 14
    end

end
