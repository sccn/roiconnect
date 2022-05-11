function [out] = DTF(A,Sigma,f_array,type)
% DTF - Directed Transfer Function (formula (11) from Baccal? and Sameshima
% 2001)
%   [Y] = DTF(A,C,TYPE) computes the DTF given the MVAR coefficients matrix
%   A and the residuals covariance matrix C. TYPE determines whether to
%   compute the standard DTF (TYPE='DTF') or the Directed DTF (TYPE = 'dDTF')
%   or the ffDTF (TYPE = 'ffDTF')
%
%   See also: PDC

%   Author: G. Gomez-Herrero, german.gomezherrero@ieee.org
%   $Revision: 1.0 $Date: 22/06/2007

if nargin < 4, type = 'dtf'; end


[N,L] = size(A);
p = L/N;
A = reshape(A,[N,N,p]);
densum = zeros(N,1);

% initialize variables
if strcmpi(type,'dtf'),
    dtf = zeros(N,N,length(f_array));
elseif sum(strcmpi(type,{'ffdtf','ddtf'})),
    ffdtf = zeros(N,N,length(f_array));
    S = zeros(N,N,length(f_array));
else
    error('Wrong type');
end
for i = 1:length(f_array)
    f = f_array(i);
    % build A(f)
    Af = zeros(N,N);
    for r = 1:p
        Af = Af + squeeze(A(:,:,r))*exp(-j*2*pi*f*r);
    end
    Hf = inv(eye(N,N)-Af);

    den = sum(((abs(Hf).^2)),2);
    if strcmpi(type,'dtf'),
        dtf(:,:,i) = Hf./repmat(sqrt(den),1,N);
    elseif strcmpi(type,'ffdtf') || strcmpi(type,'ddtf'),
        ffdtf(:,:,i) = Hf;
        densum = densum + den;
    end
    if strcmpi(type,'ddtf'),
        % build the power spectra matrix
        S(:,:,i) = Hf*Sigma*ctranspose(Hf);
    end
end
if strcmpi(type,'dtf'),
    out = dtf;
elseif strcmpi(type,'ffdtf'),
    out = ffdtf./repmat(sqrt(densum),[1,N,length(f_array)]);
elseif strcmpi(type,'ddtf'),
    ffdtf=ffdtf./repmat(sqrt(densum),[1,N,length(f_array)]);
    y = zeros(N,N,length(f_array));
    for i = 1:N
        for l = 1:N
            for k = 1:length(f_array)
                if abs((minor(S(:,:,k),i,i)*minor(S(:,:,k),l,l))-(minor(S(:,:,k),i,l)^2)) < eps,
                    y(i,l,k) = 1;
                else
                    y(i,l,k) = (minor(S(:,:,k),i,l)^2)/(minor(S(:,:,k),i,i)*minor(S(:,:,k),l,l));
                end
                y(i,l,k) = sqrt(y(i,l,k))*ffdtf(i,l,k);
            end
        end
    end
    out = y;
end

