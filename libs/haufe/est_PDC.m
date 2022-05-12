function [pdc] = PDC(A,f_array)
% PDC - Partial Directed Coherence Factor (formula (18) from Baccal? and
% Sameshima 2001)
%   [Y] = PDC(A,F) computes the PDC given the MVAR coefficients matrix A
%
%   See also: DTF

%   Author: G. Gomez-Herrero, german.gomezherrero@ieee.org
%   $Revision: 1.0 $Date: 21/06/2007


[N,L] = size(A);
p = L/N;
A = reshape(A,[N,N,p]);
pdc = zeros(N,N,length(f_array));
for i = 1:length(f_array)
    f = f_array(i);
    % build A(f)
    Af = zeros(N,N);
    for r = 1:p
        Af = Af + squeeze(A(:,:,r))*exp(-j*2*pi*f*r);
    end
    Afh = (eye(N,N)-Af);
    den = zeros(1,N);
    for k = 1:N,
        den(1,k) = sqrt(ctranspose(Afh(:,k))*Afh(:,k));
    end
    pdc(:,:,i) = Afh./repmat(den,N,1);
end


