%% FILTER A VECTOR NOISE WITH A SPECIFIED STRICTLY CAUSAL MVAR MODEL: Y(n)=A(1)Y(n-1)+...+A(p)Y(n-p)+U(n)

%%% INPUT
% A=[A(1)...A(p)]: M*pM matrix of the MVAR model coefficients (strictly causal model)
% U: M*N matrix of innovations

%%% OUTPUT
% Y: M*N matrix of simulated time series

function [Y]=MVARfilter(A,U)

N=length(U);
M=size(A,1);
p=size(A,2)/M;

% Y(n)=A(1)Y(n-1)+...+A(p)Y(n-p)+U(n)
Y=zeros(M,N);
for n=1:N
    for k=1:p
        if n-k<=0, break; end; % if n<=p, stop when k>=n
        Y(:,n)=Y(:,n) + ( A(:,(k-1)*M+(1:M)) * Y(:,n-k) );
    end
    Y(:,n)=Y(:,n)+U(:,n);
end
