%% realization of the instantaneous model : U = L*W
%%% OUTPUT
% U: N*M matrix of filtered noises

% INPUT
% N data length
% C: input covariance matrix (may be interpreted as Su or Sw, see above)
% B0: M*M matrix of instantaneous effects (when relevant)

% when flag='StrictlyCausal':
% given Su, applies Cholesky decomposition to find L and Sw
% then generates U = L*W, for a realization of gaussian W of variance Sw

% when flag='ExtendedGauss':
% given Sw and B(0), computes L=[I-B(0)]^(-1)
% then generates U = L*W, for a realization of gaussian W of variance Sw

% when flag='ExtendedNonGauss':
% given Swand B(0), computes L=[I-B(0)]^(-1)
% then generates U = L*W, for a realization of nongaussian W of variance Sw


function U=InstModelfilter(N,C,flag,B0)

error(nargchk(3,4,nargin));%min and max input arguments


M=size(C,1);


switch flag
case {'StrictlyCausal'} % C is Su
    [L,Sw]=choldiag(C);
    W = randn(M,N); % W independent and gaussian
    for m=1:M % This normalizes W to have the appropriate variance (and zero mean)
        W(m,:)=sqrt(Sw(m,m))*(W(m,:)-mean(W(m,:)))/std(W(m,:));
    end
    U=L*W;
        
case {'ExtendedGauss'} % C is Sw
    invL=eye(M)-B0;
    if det(invL)==0, error('B0 is not invertible, ill-conditioned problem!'), end;
    L=inv(invL);
    W = randn(M,N); % W independent and gaussian
    for m=1:M % This normalizes W to have the appropriate variance (and zero mean)
        W(m,:)=sqrt(C(m,m))*(W(m,:)-mean(W(m,:)))/std(W(m,:));
    end
    U=L*W;
        
case {'ExtendedNonGauss'} % C is Sw
    invL=eye(M)-B0;
    if det(invL)==0, error('B0 is not invertible, ill-conditioned problem!'), end;
    L=inv(invL);
    %note: here we generate W independent but non-Gaussian
    % Nonlinearity exponent, selected to lie in [0.5, 0.8] or [1.2, 2.0]. (<1 gives subgaussian, >1 gives supergaussian)
    q = rand(M,1)*1.1+0.5;    
    ind = find(q>0.8);           
    q(ind) = q(ind)+0.4;     
    % This generates the disturbance variables, which are mutually independent, and non-gaussian
    W = randn(M,N);
    W = sign(W).*(abs(W).^(q*ones(1,N)));
    % This normalizes the disturbance variables to have the appropriate scales
    W = W./( (  sqrt(mean((W').^2)') ./ sqrt(diag(C)) )*ones(1,N) );
    U=L*W;        
        
end
