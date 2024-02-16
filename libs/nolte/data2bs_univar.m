
function [cs,csnr,nave]=data2bs_univar(data,segleng,segshift,epleng,maxfreqbins,para)
% calculates bispectrum  from data in general for event-related measurement
% as univariate measures, i.e., always within each sensor or within a fixed sensor combination. 
%
% usage: [cs,csnr,nave]=data2bs_univar(data,segleng,segshift,epleng,maxfreqbins,para);
%
% input: 
% data: ndat times nchan matrix each colum is the time-series in one
%             channel;
% segleng: length of each segment in bins, e.g. segleng=1000;  
% segshift: numer of bins by which neighboring segments are shifted;
%           e.g. segshift=segleng/2 makes overlapping segments
% epleng: leng of each epoch
% maxfreqbins: maximum frequency in bins, starting at zeros Hertz. 
%             The frequency resolution, df, is given by the physical length of a
%             segment, say T. Then df=1/T. E.g. if T=2 seconds, the maxfreqbins=101
%             means that the maximum physical frequency is 50 Hertz.
%
% if para.chancomb = [i,j,k], the bispectrum is computed for the
%   respective channel combination. For example, if para.chancomb = [1,2,1],
%   the output will be cs([1,2,1], f1, f2) = <x(f1)_1*x(f2)_2*conj(x(f1+f2-1)_1)>
%   (for a comparison, check the default output). In this example, at least two channels 
%   are required.
%
% output: 
% cs: nchan  by nf by nf tensor for nf frequencies (i.e. nf=maxfreqbins)   
%  cs(i,f1,f2)=<x(f1)_i*x(f2)_i*conj(x(f1+f2-1)_i)>
%  where  x is the Fourier-transform of the data of each segment
%
% csn:  corresponding normalization factor defined by 
%       csn(i,f1,f2)=N_i(f1) N_i(f2) N_i(f1+f2-1);
%       where N_p(f) is defined as (<abs(x(f)_p)^3>)^(1/3) 
%       Bicoherence can be calculated as cs./csn 
%
% nave: number of averages

[ndat,nchan]=size(data);
nf=maxfreqbins;


mywindow=repmat(hanning(segleng),1,nchan);
if nargin>5
    if isfield(para,'fave')
     % fave=para.fave;
     end;
    if isfield(para,'mywindow');
       mywindow=repmat(para.mywindow,1,nchan);
    end
    if isfield(para,'chancomb')
        try
            ichan = para.chancomb(1);
            jchan = para.chancomb(2);
            kchan = para.chancomb(3);

            % check if there are enough channels
            max_chans = max([ichan, jchan, kchan]);
            if nchan < max_chans
                warning('Not enough channels, check the documentation of the "para.chancomb" parameter. Proceeding with the default option.')
                para = rmfield(para, 'chancomb');
            end
        catch
            warning('Problem using the "para.chancomb" parameter. Please check the documentation. Proceeding with the default option.')
            para = rmfield(para, 'chancomb');
        end
    end
end

nep=floor(ndat/epleng);
nseg=floor((epleng-segleng)/segshift)+1; %total number of segments

if nargin > 5 && isfield(para, 'chancomb')
    cs = zeros(nf, nf);
    csnr = zeros(nf, nf);
    csn = zeros(1, 2*nf-1);
else
    cs=zeros(nchan,nf,nf);
    csnr=zeros(nchan,nf,nf);
    csn=zeros(nchan,2*nf-1);
    %csloc=zeros(nf,nf);
end

%figure;plot(mywindow);
nave=0;
for j=1:nep;
    %disp(j)
    dataep=data((j-1)*epleng+1:j*epleng,:);
    for i=1:nseg; %average over all segments;
      dataloc=dataep((i-1)*segshift+1:(i-1)*segshift+segleng,:);
      datalocfft=fft(dataloc.*mywindow);
      datalocfft=datalocfft(1:2*nf-1,:);
      cslocn=((abs(datalocfft)).^3)';
      if nargin > 5 && isfield(para, 'chancomb')
          xx=hankel(conj(datalocfft(1:2*nf-1,kchan)));
          csloc(:,:)=(datalocfft(1:nf,ichan)*transpose(datalocfft(1:nf,jchan))).*xx(1:nf,1:nf);
          cs(:,:)=cs(:,:)+csloc(:,:);
      else
          for ichan=1:nchan;
              xx=hankel(conj(datalocfft(1:2*nf-1,ichan)));
              csloc(ichan,:,:)=(datalocfft(1:nf,ichan)*transpose(datalocfft(1:nf,ichan))).*xx(1:nf,1:nf);
              cs(ichan,:,:)=cs(ichan,:,:)+csloc(ichan,:,:);
          end
      end
      %xxx=squeeze(csloc(1,1:3,1:3))
      
      nave=nave+1;
      csn=csn+cslocn;
    end;
end

cs=cs/nave;
csn=csn/nave;
csn=power(csn,1/3);
for i=1:nchan;for f1=1:nf;for f2=1:nf;
    csnr(i,f1,f2)=(csn(i,f1)*csn(i,f2)*csn(i,f1+f2-1));
end;end;end


  
    

return;
    
