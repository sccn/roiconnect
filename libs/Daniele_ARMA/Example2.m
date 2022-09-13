%% Example2 Stokes PNAS2017 - computation of DC
clear; close all; clc;

%% Parameters
fres=1000;  % frequency resolution for spectral causalities
jj=2; % target
ii=1; % driver

%%%% Simulation parameters
fc=120; % sampling freq in Hz
M=2;
plag=2;
r1=0.65; f1o=50; f1=f1o/fc; % f1o=10, 20, 30, 40, 50 Hz
r2=0.65; f2o=[10 30 50];% f2o=10,30,50 Hz
Su=eye(M);


%% analysis cycle
for is=1:length(f2o)
    f2=f2o(is)/fc;
    % theoretical coeffs
    Ak=zeros(M,M,plag);
    for m=1:M
        eval(['Ak(m,m,1)=2*r' int2str(m) '*cos(2*pi*f' int2str(m) ');']);
        eval(['Ak(m,m,2)=-r' int2str(m) '^2;']);
    end
    Ak(2,1,1)=1;
    Am=[];
    for kk=1:plag
        Am=[Am Ak(:,:,kk)];
    end
    
    %%% Spectral functions evaluated using the blockMVAR toolbox
    Mv=[1 1];
    [bDC,bPDC,mF,mG,bS,bP,bH,bAf,freq] = block_fdMVAR(Am,Su,Mv,fres,fc);
    f(:,is)=squeeze(mF(jj,ii,:)); % in this bivariate case, this is equal to the Granger-Geweke spectral causality 
    Sjj(:,is)=abs(cell2mat(squeeze(bS(jj,jj,:)))); % spectrum of receiver
    Sii(:,is)=abs(cell2mat(squeeze(bS(ii,ii,:)))); % spectrum of transmitter
    DCji(:,is)=squeeze(bDC(jj,ii,:)); % directed coherence transmitter->receiver
    DCjj(:,is)=squeeze(bDC(jj,jj,:)); % directed coherence receiver->receiver
    
end
% Spectral decomposition
Sj_i=DCji.*Sjj;
Sj_j=DCjj.*Sjj;
SjjdB=10*log10(Sjj);


%% Figure
set(0,'defaultAxesFontSize', 13);
set(0,'defaultTextFontSize', 14);
colpS=[146 54 135]/255;
coldc=[64 128 0]/255;
colnondc=[255 128 0]/255;

figure(2); clf
for is=1:3
    subplot(3,4,(is-1)*4+1);
    plot(freq,Sii(:,is),'Color',coldc,'LineWidth',1.5);
    xlim([0 fc/2]); 
    set(gca,'XTick',[0,fc/4,fc/2]);
    set(gca,'XTickLabel',{'0',num2str(fc/4),num2str(fc/2)});
    ylabel(['S_{' int2str(ii) int2str(ii) '}'],'FontSize',16);
    xlabel('Frequency (Hz)'); set(gca,'XTick',[0,10,30,50,60]);
    set(gca,'XTickLabel',{'0','10','30','50','60'});
    if is==1,
         title('Spectrum transmitter');
    end
    
    subplot(3,4,(is-1)*4+2);
    plot(freq,Sjj(:,is),'Color',colpS,'LineWidth',1.5);
    hold on;
    plot(freq,Sj_j(:,is),'Color',colnondc,'LineWidth',1.5);
    plot(freq,Sj_i(:,is),'Color',coldc,'LineWidth',1.5); % superimpose the version of GG causality
    xlim([0 fc/2]);
    set(gca,'XTick',[0,fc/4,fc/2]);
    set(gca,'XTickLabel',{'0',num2str(fc/4),num2str(fc/2)});
    ylabel(['S_{' int2str(jj) int2str(jj) '}'],'FontSize',16);
    if is==1, title('Spectrum receiver');
    hleg=legend(['S_{' int2str(jj) int2str(jj) '}(f)'], ['S_{' int2str(jj) '|' int2str(jj) '}(f)'],['S_{' int2str(jj) '|' int2str(ii) '}(f)']);
    set(hleg,'FontSize',12); %set(hleg,'visible','off');
    legend boxoff;
    end
    xlabel('Frequency (Hz)'); set(gca,'XTick',[0,10,30,50,60]);
    set(gca,'XTickLabel',{'0','10','30','50','60'});
    
    subplot(3,4,(is-1)*4+3);
    plot(freq,DCji(:,is),'Color',coldc,'LineWidth',1.5);
    xlim([0 fc/2]);
    set(gca,'XTick',[0,fc/4,fc/2]);
    set(gca,'XTickLabel',{'0',num2str(fc/4),num2str(fc/2)});
    ylabel(['DC_{' int2str(ii) '\rightarrow' int2str(jj) '}'],'FontSize',16);
    if is==1,
        title('directed coherence');
    end
    xlabel('Frequency (Hz)'); set(gca,'XTick',[0,10,30,50,60]);
    set(gca,'XTickLabel',{'0','10','30','50','60'});
    
    subplot(3,4,(is-1)*4+4);
    plot(freq,f,'Color',colpS,'LineWidth',1.5);
    xlim([0 fc/2]);
    set(gca,'XTick',[0,fc/4,fc/2]);
    set(gca,'XTickLabel',{'0',num2str(fc/4),num2str(fc/2)});
    ylabel(['GG_{' int2str(ii) '\rightarrow' int2str(jj) '}'],'FontSize',16);
    if is==1,
        title('spectral causality');
    end
    xlabel('Frequency (Hz)'); set(gca,'XTick',[0,10,30,50,60]);
    set(gca,'XTickLabel',{'0','10','30','50','60'});
    
end

%%%% use export_fig toolbox to export
% set(gcf, 'Position', [1 1 1400 900]);
% set(gcf, 'Color', 'w');
% export_fig test2.png -m2.5






