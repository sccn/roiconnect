%% CHOLESKY DECOMPOSITION WITH DIAGONAL MATRIX
% given A  positive definite symmetric matrix, returns a lower triangular
% matrix L and a diagonal matrix D such that A=LDL'

function [L,D]=choldiag(A);

error(nargchk(1,1,nargin));%min e max di input arguments

n=size(A,1);
if sum(sum(A'~=A)), error('A non è simmetrica'); end
if size(A,1)~=size(A,2), error('A non è quadrata'); end
chol(A); %contains an implicit error implicito if A is not positive definite

L=eye(n,n);
d=zeros(n,1);
for i=1:n
   
    for j=1:i-1
        tmp=0;
        for k=1:j-1
            tmp=tmp+L(i,k)*L(j,k)*d(k);
        end
        L(i,j)=(1/d(j))*(A(i,j)-tmp);
    end
    
    tmp2=0;
    for k=1:i-1
        tmp2=tmp2+L(i,k)^2*d(k);
    end
    d(i)=A(i,i)-tmp2;
    
end
D=diag(d);

% verifica
% % L
% % D
% % L*D*L'     

        




