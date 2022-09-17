%% IDENTIFICATION OF STRICTLY CAUSAL MVAR MODEL: Y(n)=A(1)Y(n-1)+...+A(p)Y(n-p)+U(n)
% makes use of autocovariance method (vector least squares)

%%% input:
% Y, M*N matrix of time series (each time series is in a row)
% p, model order
% Mode, determines estimation algorithm (0:builtin least squares, else other methods [see mvar.m from biosig package])

%%% output:
% Am=[A(1)...A(p)], M*pM matrix of the estimated MVAR model coefficients
% S, estimated M*M input covariance matrix
% Yp, estimated time series
% Up, estimated residuals
% Z, observation matrix (often optional, useful e.g. for resampling)


function [Am,S,Yp,Up,Z,Yb]=idMVAR(Y,p,Mode)

% error(nargchk(1,3,nargin));
% if nargin < 3, Mode=0; end % default use least squares estimate
% if nargin < 2, p=10; end % default model order

[M,N]=size(Y);

%% IDENTIFICATION
Z=NaN*ones(p*M,N-p); % observation matrix
for j=1:p
    for i=1:M
        Z((j-1)*M+i,1:N-p)=Y(i, p+1-j:N-j);
    end
end

if Mode==0
    Yb=NaN*ones(M,N-p); % Ybar
    for i=1:M
        Yb(i,1:N-p)=Y(i,p+1:N);
    end
    Am=Yb/Z; % least squares!
%     fprintf('using least squares\n');        
else
    Am = mvar(Y', p, Mode); % estimates from biosig code
%     fprintf(['using biosig ' int2str(Mode) ' mode\n']);  
end

Yp=Am*Z; 
Yp=[NaN*ones(M,p) Yp]; % Vector of predicted data

Up=Y-Yp; Up=Up(:,p+1:N); % residuals of strictly causal model
S=cov(Up');


