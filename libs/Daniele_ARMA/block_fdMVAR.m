%% FREQUENCY DOMAIN BLOCK MVAR ANALYSIS
% References:
% L.Faes and G. Nollo, "Measuring Frequency Domain Granger Causality for Multiple Blocks of Interacting Time Series", Biological Cybernetics 2013. DOI: 10.1007/s00422-013-0547-5
% L.Faes, S. Erla and G. Nollo, "Block Partial Directed Coherence: a New Tool for the Structural Analysis of Brain Networks", International Journal of Bioelectromagnetism, Vol. 14, No. 4, pp. 162 - 166, 2012

%%% inputs: 
% Am=[A(1)...A(p)]: Q*pQ matrix of the MVAR model coefficients (strictly causal model)
% Su: Q*Q covariance matrix of the input noises
% Mv: number of series in each block
% N= number of points for calculation of the spectral functions (nfft)
% Fs= sampling frequency

%%% outputs:
% bDC= block Directed Coherence (Eq. 14a)
% bPDC= block Partial Directed Coherence (Eq. 14b)
% mF= multivariate total causality (Eq. 13a)
% mG= multivariate direct causality (Eq. 13b)
% bS= block spectral density matrix (Eq. 12a)
% bP= inverse block spectral density matrix (Eq. 12b)
% bH= block transfer matrix
% bAf= block spectral coefficient matrix
% f= vector of frequencies


function [bDC,bPDC,mF,mG,bS,bP,bH,bAf,f] = block_fdMVAR(Am,Su,Mv,N,Fs)
% clear; close all; clc;
% [Bm,B0,Sw]=simuMVARcoeff(1);
% Am=Bm;Su=Sw;
% Mv=[2 1 1]'; % vector of Mi (dimension of each block)
% N=512;
% Fs=1;
% f = (0:N-1)*(Fs/(2*N));

%%
Q= size(Am,1); % Am has dim Q*pQ
M=length(Mv);
p = size(Am,2)/Q; % p is the order of the MVAR model

if nargin<5, Fs= 1; end;   
if nargin<4, N = 512; end;
if all(size(N)==1),	 %if N is scalar
    f = (0:N-1)*(Fs/(2*N)); % frequency axis
else            % if N is a vector, we assume that it is the vector of the frequencies
    f = N; N = length(N);
end;

s = exp(sqrt(-1)*2*pi*f/Fs); % vector of complex exponentials
z = sqrt(-1)*2*pi/Fs;


%% Initializations: spectral matrices have M rows, M columns and are calculated at each of the N frequencies
A = [eye(Q) -Am]; % matrix from which M*M blocks are selected to calculate spectral functions
invSu=inv(Su);
% i-j block of Su and Su^(-1)
bSu=cell(M,M); binvSu=cell(M,M);  
for i=1:M
    for j=1:M
        i1=sum(Mv(1:i)); i0=i1-Mv(i)+1;
        j1=sum(Mv(1:j)); j0=j1-Mv(j)+1;
        bSu{i,j}=Su(i0:i1,j0:j1);
        binvSu{i,j}=invSu(i0:i1,j0:j1);
    end
end

% whole QxQ matrices
Af=zeros(Q,Q,N); % Coefficient Matrix in the frequency domain (it is Abar(w))
H=zeros(Q,Q,N); % Transfer Matrix H(w)
S=zeros(Q,Q,N); % Spectral Matrix S(w)
P=zeros(Q,Q,N); % Inverse Spectral Matrix P(w)
% corresponding cell array of matrices (blocks of size Mi x Mj)
bAf=cell(M,M,N);
bH=cell(M,M,N); 
bS=cell(M,M,N); 
bP=cell(M,M,N);
% causality functions
bDC=zeros(M,M,N); % block directed coherence
bPDC=zeros(M,M,N); % block partial directed coherence
mF=NaN*ones(M,M,N); % multivariate logarithmic causality f
mG=NaN*ones(M,M,N); % multivariate logarithmic direct causality g

%%% computation of spectral functions
for n=1:N, % at each frequency  
        %%% Coefficient matrix in the frequency domain
        As = zeros(Q,Q); % matrix As(z)=I-sum(A(k))
        for k = 1:p+1,
            As = As + A(:,k*Q+(1-Q:0))*exp(-z*(k-1)*f(n));  %indicization (:,k*Q+(1-M:0)) extracts the k-th Q*Q block from the matrix B (A(1) is in the second block, and so on)
        end;        
        Af(:,:,n) = As; %%% Coefficient matrix
        H(:,:,n)  = inv(As); %%% Transfer matrix     
        S(:,:,n)  = H(:,:,n)*Su*H(:,:,n)'; %%% Spectral matrix -   ' stands for Hermitian transpose
        P(:,:,n) = inv(S(:,:,n)); %%% Inverse Spectral matrix P(:,:,n) = As'*invSu*As;
        
        %%% extraction of blocks indexes
        for i=1:M
            for j=1:M
                %indexes of the i-j block
                i1=sum(Mv(1:i)); i0=i1-Mv(i)+1;
                j1=sum(Mv(1:j)); j0=j1-Mv(j)+1;
                % i-j block of all matrices of interest
                bAf{i,j,n}=Af(i0:i1,j0:j1,n);
                bH{i,j,n}=H(i0:i1,j0:j1,n);
                bS{i,j,n}=S(i0:i1,j0:j1,n);
                bP{i,j,n}=P(i0:i1,j0:j1,n);               
            end
        end
        
        % computation of causality measures at frequency n
        for i=1:M
            for j=1:M
                                    
%                 if det(bP{j,j,n}) < -0.000001, error('determinante negativo!'); end
%                 if det(bS{i,i,n}) < -0.000001, error('determinante negativo!'); end
%                 if det(bP{j,j,n} - bAf{i,j,n}'*binvSu{i,i}*bAf{i,j,n}) < -0.000001, error('determinante negativo!'); end
%                 if det(bS{i,i,n} - bH{i,j,n}*bSu{j,j}*bH{i,j,n}') < -0.000001, error('determinante negativo!'); end
                
                bDC(i,j,n) = 1 - abs(det(bS{i,i,n} - bH{i,j,n}*bSu{j,j}*bH{i,j,n}')) / abs(det(bS{i,i,n}));
                bPDC(i,j,n) = 1 - abs(det(bP{j,j,n} - bAf{i,j,n}'*binvSu{i,i}*bAf{i,j,n})) / abs(det(bP{j,j,n}));
                if i~=j
                    mF(i,j,n) = log( abs(det(bS{i,i,n})) / abs(det(bS{i,i,n} - bH{i,j,n}*bSu{j,j}*bH{i,j,n}')) );
                    mG(i,j,n) = log( abs(det(bP{j,j,n})) / abs(det(bP{j,j,n} - bAf{i,j,n}'*binvSu{i,i}*bAf{i,j,n})) );
                end
            end
        end      
        
end;
