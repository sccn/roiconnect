%% VARMA with B0 term to (Innovations form) State Space parameters
% computes innovations form parameters for a state space model from VARMA
% parameters using Aoki's method - this version allows for zero-lag MA coefficients

function [A,C,K,R,lambda0] = varma2iss(Am,Bm,V,B0)

%   INPUT: VARMA parameters Am, Bm, V=cov(U)
%   OUTPUT: innovations form SS parameters A, C, K, R

%%%%% internal test
%variables to be passed are Am, Bm, B0, V=Su
% clear; close all; clc;
% Am=[0.9 0 0 0.5; 0 0.6 0.2 0];
% Bm=[0.5 0; 0 0.5]; B0=Bm./5;
% V=eye(2);
%  
M = size(Am,1); %dimension of observed process
p=floor(size(Am,2)/M); %number of AR lags
q=floor(size(Bm,2)/M); %number of MA lags

L=M*(p+q); % dimension of state process (SS order)

C=[Am Bm];
R=B0*V*B0';

Ip=eye(M*p);
Iq=eye(M*q);
A11=[Am;Ip(1:end-M,:)];

if q==0
    A=A11;
    K=[eye(M); zeros(M*(p-1),M)];
else
    A12=[Bm;zeros(M*(p-1),M*q)];
    A21=zeros(M*q,M*p);
    A22=[zeros(M,M*q); Iq(1:end-M,:)];
    A=[A11 A12; A21 A22];

    K=[eye(M); zeros(M*(p-1),M); inv(B0); zeros(M*(q-1),M)];
end


% determine the variance of the process lambda0=E[Yn Yn']
O=dlyap(A,K*R*K');
lambda0=C*O*C'+R;

end


