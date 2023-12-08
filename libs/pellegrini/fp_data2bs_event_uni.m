function [cs,nave]=fp_data2bs_event_uni(data,segleng,segshift,epleng,freqpairs,nshuf)
% Calculates bispectral-tensors  from data for event-related measurement
% and their null distributions using a shuffling approach. 
%
% Based on the data2bs_event function of the METH toolbox by Guido Nolte. 
% Modified by (c) 2023 Franziska Pellegrini and Stefan Haufe
% 
% usage: [cs,nave]=data2bs_event(data,segleng,segshift,epleng,freqpairs,para);
%
% input:
% data: ndat times nchan matrix each colum is the time-series in one
%             channel;
% segleng: length of each segment in bins, e.g. segleng=1000;
% segshift: numer of bins by which neighboring segments are shifted;
%           e.g. segshift=segleng/2 makes overlapping segments
% epleng: leng of each epoch
% freqpairs: pairs of  frequency in bins
% nshuf: required number of samples (shuffles) in the null distribution 
% para: structure which is eventually used later
%
% output:
% cs: nchan by nchan by nchan by peak_combination (=2) by number of shuffles 
% tensor such that
%  cs(i,j,k,ishuf)=<x(f1)_i*x(f2)_j*conj(x(f1+f2-1)_k)>
%  where f1=freqpairs(f,1) and  f2=freqpairs(f,2),
%  x=fft(data) and the average is over epeochs and segments. 
%
% NOTE that this function is restricted to one frequency combination only!
% The addition for frequency BANDS will be added in the future.
%
% if para.fave=0 then cs contains a fifth argument denoting
% the   segment.
% if para.fave=1 or ommited, then cs was averaged over segments.

% nave: number of averages

[ndat,nchan]=size(data);

nep=floor(ndat/epleng);

nseg=floor((epleng-segleng)/segshift)+1; %total number of segments
assert(nseg==1,'only possible with 1 segment')

cs=zeros(nchan,nchan,nchan,2,nshuf);

%get Fourier coefficients
coeffs = fp_fft_coeffs(data,segleng,segshift,epleng,freqpairs);

for ishuf = 1:nshuf
    nave=0;
    csloc1=zeros(nchan,nchan,nchan);
    csloc2=zeros(nchan,nchan,nchan);
    cs1=zeros(nchan,nchan,nchan);
    cs2=zeros(nchan,nchan,nchan);

    if ishuf == 1 %the first shuffle contains the true order
        inds = 1:nep;
    else
        inds = randperm(nep,nep); %indices for shuffling of epochs for channel 2 
    end
    
    for j=1:nep %loop over epochs 
        
        %bispec of f1 (low peak), f2-f1 (left side lobe), f2 (high peak)
        csloc1(1,1,1)=transpose(coeffs(1,1,j))        *coeffs(2,1,j)      *conj(coeffs(3,1,j));
        csloc1(2,1,1)=transpose(coeffs(1,2,inds(j)))  *coeffs(2,1,j)      *conj(coeffs(3,1,j));
        csloc1(1,2,1)=transpose(coeffs(1,1,j))        *coeffs(2,2,inds(j))*conj(coeffs(3,1,j));
        csloc1(1,1,2)=transpose(coeffs(1,1,j))        *coeffs(2,1,j)      *conj(coeffs(3,2,inds(j)));
        csloc1(2,2,1)=transpose(coeffs(1,2,inds(j)))  *coeffs(2,2,inds(j))*conj(coeffs(3,1,j));
        csloc1(1,2,2)=transpose(coeffs(1,1,j))        *coeffs(2,2,inds(j))*conj(coeffs(3,2,inds(j)));
        csloc1(2,2,2)=transpose(coeffs(1,2,inds(j)))  *coeffs(2,2,inds(j))*conj(coeffs(3,2,inds(j)));
        csloc1(2,1,2)=transpose(coeffs(1,2,inds(j)))  *coeffs(2,1,j)      *conj(coeffs(3,2,inds(j)));
        
        %bispec of f1 (low peak), f2 high peak), f1+f2 (right side lobe)
        csloc2(1,1,1)=transpose(coeffs(1,1,j))          *coeffs(3,1,j)      *conj(coeffs(4,1,j));
        csloc2(2,1,1)=transpose(coeffs(1,2,inds(j)))    *coeffs(3,1,j)      *conj(coeffs(4,1,j));
        csloc2(1,2,1)=transpose(coeffs(1,1,j))          *coeffs(3,2,inds(j))*conj(coeffs(4,1,j));
        csloc2(1,1,2)=transpose(coeffs(1,1,j))          *coeffs(3,1,j)      *conj(coeffs(4,2,inds(j)));
        csloc2(2,2,1)=transpose(coeffs(1,2,inds(j)))    *coeffs(3,2,inds(j))*conj(coeffs(4,1,j));
        csloc2(1,2,2)=transpose(coeffs(1,1,j))          *coeffs(3,2,inds(j))*conj(coeffs(4,2,inds(j)));
        csloc2(2,2,2)=transpose(coeffs(1,2,inds(j)))    *coeffs(3,2,inds(j))*conj(coeffs(4,2,inds(j)));
        csloc2(2,1,2)=transpose(coeffs(1,2,inds(j)))    *coeffs(3,1,j)      *conj(coeffs(4,2,inds(j)));
        
        cs1=cs1+csloc1;
        cs2=cs2+csloc2;
        
        nave=nave+1;
    end
    
    %shape cs: chan x chan x chan x peak_combination x shuffles 
    cs(:,:,:,1,ishuf) = cs1./nave; 
    cs(:,:,:,2,ishuf) = cs2./nave; 
    
end