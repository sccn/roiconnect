function [K,V,rep,L,P] = ss2iss(A,C,Q,R,S)

% Compute innovations form parameters for a state space model in general form by
% solution of a discrete algebraic Riccati equation (DARE) (eqs. 7, 8a, 8b in
% the reference article).
%
% A,C,Q,R,S - general form state space parameters
%
% K         - Kalman gain matrix
% V         - innovations covariance matrix
% rep       - DARE report (see below)
% L         - DARE stablising eigenvalues
% P         - DARE solution
%
% The value returned in rep is negative if an unrecoverable error was detected:
% rep = -1 means that the DARE was found to have eigenvalues on (or near) the
% unit circle, while rep = -2 indicates that no stabilising solution to the DARE
% could be found. If no error occurred, rep returns the relative residual of the
% DARE solution, which should be tested (rep > sqrt(eps) is a reasonable test
% for accuracy problems).

[m,m1]  = size(A); assert(m1 == m);
[n,m1]  = size(C); assert(m1 == m);
[m1,m2] = size(Q); assert(m1 == m && m2 == m);
[m1,n1] = size(S); assert(m1 == m && n1 == n);
[n1,n2] = size(R); assert(n1 == n && n2 == n);

K = [];
V = [];
P = [];

% We solve the DARE using the Generalized Schur (QZ) decomposition method

% Extended pencil

H = [A' zeros(m) C'; -Q  eye(m) -S; S' zeros(n,m) R];
J = [eye(m) zeros(m,m+n); zeros(m) A zeros(m,n); zeros(n,m) -C zeros(n)];

% NOTE - we don't balance!

mm = 2*m;
[q,~] = qr(H(:,mm+1:mm+n));
H = q(:,n+1:mm+n)'*H(:,1:mm);
J = q(:,n+1:mm+n)'*J(:,1:mm);

% QZ algorithm

realHJ = isreal(H) && isreal(J);
i = 1:mm;
if realHJ
    [JJ,HH,q,z] = qz(J(i,i),H(i,i),'real');
else
    [JJ,HH,q,z] = qz(J(i,i),H(i,i),'complex');
end
[JJ,HH,~,z(i,:),qzflag] = ordqz(JJ,HH,q,z,'udo');
L = ordeig(JJ,HH);

% Check for stable invariant subspace

sis = abs(L) > 1;
if ~qzflag || any(~sis(1:m,:)) || any(sis(m+1:mm,:))
    rep = -1;
    return % IMPORTANT: caller must test!!! (error is: ''DARE: eigenvalues on/near unit circle'')
end

P1 = z(1:m,1:m);
P2 = z(m+1:mm,1:m);

% Solve P*P1 = P2

[LL,UU,pvec] = lu(P1,'vector');
if rcond(UU) < eps
    rep = -2;
    return % IMPORTANT: caller must test!!! (error is: 'DARE: couldn''t find stablising solution')
end
P(:,pvec) = (P2/UU)/LL;
P = (P+P')/2;

% Compute Kalman gain matrix K and innovations covariance matrix V

U = A*P*C'+S;
V = C*P*C'+R;
K = U/V;

% Check accuracy

APA = A*P*A'-P;
UK = U*K';
rep = norm(APA-UK+Q,1)/(1+norm(APA,1)+norm(UK,1)+norm(Q,1)); % relative residual

% IMPORTANT: test for accuracy  - something like
%
% if rep > sqrt(eps)
%     warning('DARE: there were accuracy issues (relative residual = %e)',rep);
% end

if nargout > 3

 % Return stable eigenvalues

    if realHJ
        L = L(m+1:mm);
    else
        L = conj(L(m+1:mm));
    end
end
