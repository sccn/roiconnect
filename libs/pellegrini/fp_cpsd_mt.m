function S = fp_cpsd_mt(X1,X2,ind_1,ind_2,h,window,noverlap,nchunks,taparray)

% Copyright (c) 2022 Franziska Pellegrini and Stefan Haufe

ind_pow = intersect(ind_1, ind_2);
nfft = 2*(h-1);

n1 = numel(ind_1);
n2 = numel(ind_2);
n3 = numel(ind_pow);
S = complex(zeros(h,n1,n2));
pow = zeros(h,n3);

winstep = window-noverlap;
ntapers = size(taparray,2);

% compute tapered periodogram with FFT

for k = 1:nchunks
    
    XSEG1 = X1((1:window) + (k-1)*winstep,:);
    XSEG2 = X2((1:window) + (k-1)*winstep,:);
    
    % compute periodogram
    P1 = fft(taparray.*permute(XSEG1(:,:,ones(1,ntapers)),[1 3 2]),nfft);
    P1 = P1(1:h,:,:);
    P2 = fft(taparray.*permute(XSEG2(:,:,ones(1,ntapers)),[1 3 2]),nfft);
    P2 = P2(1:h,:,:);
    
    % now make cross-products of them to fill cross-spectrum matrix
    
    for ii = 1:n1
        o = ind_1(ii); 
        for jj = 1:n2
            oo = ind_2(jj);      
            S(:,ii,jj) = S(:,ii,jj) + mean(P1(:,:,o) .* conj(P2(:,:,oo)),2);  
        end
    end
    
    if ~isempty(ind_pow)
        for pp = 1:n3 
            u = ind_pow(pp);
            pow(:,pp) = pow(:,pp) + mean(P1(:,:,u) .* conj(P1(:,:,u)),2);
        end
    end
   
end

S = S/nchunks;

if ~isempty(ind_pow)
    pow = pow/nchunks;

    for pp = 1:n3
        clear a b 
        a = find(ind_1 == ind_pow(pp));
        b = find(ind_2 == ind_pow(pp));
        S(:,a,b) = pow(:,pp);
    end
end

