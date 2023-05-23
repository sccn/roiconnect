function [norms,nave]=data2bs_threenorm(data,segleng,segshift,epleng,freqpairs,para)
% calculates bispectral-tensors  from data in general for event-related
% measurements
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
maxfreqbin=sum(freqpairs)-1;
f1=freqpairs(1);f2=freqpairs(2);
mywindow=repmat(hanning(segleng),1,nchan);
kontrand=0;
nrun=1;
if nargin>5
    if isfield(para,'fave')
        fave=para.fave;
    end
    if isfield(para,'mywindow')
        mywindow=repmat(para.mywindow,1,nchan);
    end
    
end

norms=zeros(nchan,nchan,nchan);

nep=floor(ndat/epleng);

nseg=floor((epleng-segleng)/segshift)+1; %total number of segments
datafft=zeros(maxfreqbin,nchan,nseg,nep);
norm1=zeros(nchan,1);
norm2=zeros(nchan,1);
norm3=zeros(nchan,1);
norms=zeros(nchan,nchan,nchan);

for j=1:nep
    dataep=data((j-1)*epleng+1:j*epleng,:);
    for i=1:nseg %average over all segments;
        dataloc=dataep((i-1)*segshift+1:(i-1)*segshift+segleng,:);
        datalocfft=fft(detrend(dataloc).*mywindow);
        datafft(:,:,i,j)=datalocfft(1:maxfreqbin,:);
    end
end


nave=0;
for j=1:nep
    for i=1:nseg
        for k=1:nchan
            norm1(k)=norm1(k)+abs(datafft(f1,k,i,j)).^3;
            norm2(k)=norm2(k)+abs(datafft(f2,k,i,j)).^3;
            norm3(k)=norm3(k)+abs(datafft(f1+f2-1,k,i,j)).^3;
        end
        nave=nave+1;
    end   
end

norm1=(norm1/nave).^(1/3);
norm2=(norm2/nave).^(1/3);
norm3=(norm3/nave).^(1/3);
for i1=1:nchan
    for i2=1:nchan
        for i3=1:nchan
            norms(i1,i2,i3)=norm1(i1)*norm2(i2)*norm3(i3);
        end
    end
end
