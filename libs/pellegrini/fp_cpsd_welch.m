function S = fp_cpsd_welch(X_1, X_2,ind_1,ind_2,h,window,noverlap)
%X_1 and X_2 can be the same - then the usual CS is
%calculated. Otherwise X_2 should contain data from another trial than
%X_1. ind_1 and ind_2 contain channels of interest. The output S will be of 
%size length(ind_1) x length(ind_2) x nfreq. 

% Copyright (c) 2022 Franziska Pellegrini and Stefan Haufe

ind_pow = intersect(ind_1, ind_2);
nfft = 2*(h-1);

n1 = numel(ind_1);
n2 = numel(ind_2);
S = complex(zeros(n1,n2,h));

for ii = 1:n1
    o = ind_1(ii); 
    S(ii,:,:) = transpose(cpsd(X_1(:,o),X_2(:,ind_2),window,noverlap,nfft)); % cross-spectra 
    if ismember(o,ind_pow)
        clear b 
        b = find(ind_2 == o);        
        S(ii,b,:) = transpose(cpsd(X_1(:,o),X_1(:,o),window,noverlap,nfft));
    end 
    
end       
        