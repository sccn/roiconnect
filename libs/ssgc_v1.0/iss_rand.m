function [A,C,K,V] = iss_rand(n,m,rhoa,dis)

% Generate some random stable and minimum-phase normally distributed innovations
% form parameters.
%
% n     - observation variable dimension
% m     - state dimension
% rhoa  - state transition matrix spectral norm
% dis  - display spectral norm of pencil of AR transition matrices
%
% A,C,K - innovations form state space parameters
%
% We need rhoa < 1 for a stable model. This routine creates A,K,C as random
% Gaussian matrices, and scales A so that its spectral norm is equal to rhoa. C
% and K are then each scaled by a factor sqrt(r), so that the AR transition
% matrix A - r*K*C has spectral norm < 1, thereby ensuring minimum phase. (The
% subroutine speclim is used to identify minimum and maximum r between which
% minimum phase holds, and r is chosen uniformly from this range.)

assert(rhoa < 1);

if nargin < 4 || isempty(dis), dis = false; end

A = randn(m);
A = (rhoa/specnorm(A))*A;
C = randn(n,m);
K = randn(m,n);

M = K*C;
rmin = speclim(A,M,-1,0);
rmax = speclim(A,M,+1,0);

r = rmin + (rmax-rmin)*rand; % uniform on [rmin,rmax]
sqrtr = sqrt(abs(r));
C = sqrtr*C;
K = sign(r)*sqrtr*K;

if dis
    nr = 1000;
    ramax = 1.1*max(rmax,-rmin);
    rr = linspace(-ramax,ramax,nr)';
    rrhob = zeros(nr,1);
    for i = 1:nr
        rrhob(i) = specnorm(A-rr(i)*M);
    end
    rholim = [0.9*min(rrhob) 1.1*max(rrhob)];
    plot(rr,rrhob);
    xlim([-ramax,ramax]);
    ylim(rholim);
    line([-ramax,ramax],[1 1],'Color','k');
    line([0 0],rholim,'Color','r');
    line([r r],rholim,'Color','g');
    xlabel('r');
    ylabel('rho');
    legend('rho(B)');
    title(sprintf('rho(A) = %g, rho(B) = %g',rhoa,specnorm(A-r*M)));
end

function r = speclim(A,M,r1,r2)

assert(specnorm(A-r1*M) > 1 && specnorm(A-r2*M) < 1);
while true
    r = (r1+r2)/2;
    rho = specnorm(A-r*M);
    if rho > 1, r1 = r; else r2 = r; end
    if abs(r1-r2) < eps, break; end
end

