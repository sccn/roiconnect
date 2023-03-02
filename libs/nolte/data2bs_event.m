function [cs,nave] = data2bs_event(data,segleng,segshift,epleng,freqpairs,nshuf)
% calculates bispectral-tensors  from data for event-related measurement
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
% para: structure which is eventually used later
%
% output:
% cs: nchan by nchan by nchan by number_of_frequency_pairs
%      (by number_of_segments) tensor such that
%  cs(i,j,k,f)=<x(f1)_i*x(f2)_j*conj(x(f1+f2-1)_k)>
%  where f1=freqpairs(f,1) and  f2=freqpairs(f,1),
%  x=fft(data) and the average is over epeochs and segments
%
% if para.fave=0 then cs contains a fifth argument denoting
% the   segment.
% if para.fave=1 or ommited, then cs was averaged over segments.

% nave: number of averages

[ndat,nchan]=size(data);

nep=floor(ndat/epleng);

nseg=floor((epleng-segleng)/segshift)+1; %total number of segments
assert(nseg==1,"only possible with 1 segment")

cs=zeros(nchan,nchan,nchan,2,nshuf);
nave=0;

coeffs = fp_fft_coeffs(data,segleng,segshift,epleng,freqpairs);

for ishuf = 1:nshuf
    
    csloc1=zeros(nchan,nchan,nchan);
    csloc2=zeros(nchan,nchan,nchan);
    cs1=zeros(nchan,nchan,nchan);
    cs2=zeros(nchan,nchan,nchan);

    if ishuf == 1
        inds = 1:nep;
    else
        inds = randperm(nep,nep); %indices for shuffling of epochs for f1
    end
    
    for j=1:nep
        
        for k=1:nchan
            csloc1(:,:,k)=transpose(coeffs(inds(1),:,j))*coeffs(2,:,j)*conj(coeffs(3,k,j)); %bispec of f1, f2-f1, f2
            csloc2(:,:,k)=transpose(coeffs(inds(1),:,j))*coeffs(3,:,j)*conj(coeffs(4,k,j)); %bispec of f1, f2, f1+f2
        end
        cs1=cs1+csloc1;
        cs2=cs2+csloc2;
        
        nave=nave+1;
    end
    
%     for ichan=1:nchan
%         for jchan=1:nchan
%             cs1(ichan,jchan,:)=cs1(ichan,jchan,:)/nave;
%             cs2(ichan,jchan,:)=cs2(ichan,jchan,:)/nave;
%         end
%     end
    
    cs(:,:,:,:,1,ishuf) = cs1./nave; 
    cs(:,:,:,:,2,ishuf) = cs2./nave; 
    
end