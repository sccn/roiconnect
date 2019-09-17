%% Example2 Stokes PNAS2017 - demonstrate validity of f with ISS
clear; close all; clc;

%% Parameters
numsimu=2;
N=500; %simulation length
pcrit='bic'; % 'aic', 'bic', or number for fixed model order
pmax=12; %look for BIC
lo=5; hi=95; % percentiles for showing distributions
fres=1000;  % frequency resolution for spectral causalities

jjv=[2 3 1]; % target (in Fig. 1 of PNAS)
iiv=[1 2 3]; % driver (in Fig. 1 of PNAS)

Mv=[1 1 1]; % for block analysis of conditional Geweke f measure
p1=3; p2=20; %orders according to PNAS for Geweke estimation


%%%% Simulation parameters
fc=120; % sampling freq in Hz
M=3; %generate two bivariate unidirectionally coupled AR processes,then merge
plag=3;
r1=0.9; f1o=40; f1=f1o/fc;
r2=0.7; f2o=10; f2=f2o/fc;
r3=0.8; f3o=50; f3=f3o/fc;
Su=eye(M);
% theoretical coeffs
Ak=zeros(M,M,plag);
for m=1:M
    eval(['Ak(m,m,1)=2*r' int2str(m) '*cos(2*pi*f' int2str(m) ');']);
    eval(['Ak(m,m,2)=-r' int2str(m) '^2;']);
end
Ak(2,1,1)=-0.356;
Ak(3,2,1)=-0.3098;
Ak(2,1,2)=0.7136;
Ak(3,2,2)=0.5;
Ak(2,1,3)=-0.356;
Ak(3,2,3)=-0.3098;
Am=[];
for kk=1:plag
    Am=[Am Ak(:,:,kk)];
end

% (Normalized) frequencies for frequency-domain GCs
omega = ((0:fres-1)/fres)*pi; % frequencies in range [0,pi]
z = exp(-1i*omega);         % on unit circle in complex plane
freq=fc*omega./(2*pi);


%% analysis
%%% state space model for original coeffs
Bm=[]; V=Su;
[A,C,K,R,lambda0] = varma2iss(Am,Bm,V,eye(M));

for it=1:length(jjv)
    jj=jjv(it);
    ii=iiv(it);
    kk=setdiff((1:M),[ii jj]);

    % Calculate time-domain GC
    %F = iss_GC(A,C,K,V,jj,ii);

    %%% Calculate frequency-domain GC (theoretical)
    f(:,it) = iss_SGC(A,C,K,V,z,jj,ii)';
    
    %%% Calculate frequency-domain GC (estimates)
    for is=1:numsimu
        disp(['direction ' int2str(ii) '->' int2str(jj) ', simu ' int2str(is) ' of ' int2str(numsimu)]);
        U=InstModelfilter(N,Su,'StrictlyCausal'); % gaussian innovations with covariance Su
        Y=MVARfilter(Am,U);
        
       
        % Geweke estimation (classical)
        [gewF,gewFc] = Geweke_f(Y,Mv,p1,ii,jj,fres,fc);
        ef_gew_p1(:,is,it)=gewFc;
        [gewF,gewFc] = Geweke_f(Y,Mv,p2,ii,jj,fres,fc);
        ef_gew_p2(:,is,it)=gewFc;
        
        % ISS estimation
        [eAm1,eSu1,Yp1,Up1]=idMVAR(Y,p1,0);
        [eA,eC,eK,eV,eVy] = varma2iss(eAm1,Bm,eSu1,eye(M)); % max(abs(eig(A-K*C)))
        ef(:,is,it) = iss_SGC(eA,eC,eK,eV,z,jj,ii)';
    end

end
%%% median and percentiles
ef_m=squeeze(median(ef,2));
ef_lo=squeeze(prctile(ef,lo,2));
ef_hi=squeeze(prctile(ef,hi,2));
ef_gew_p1_m=squeeze(median(ef_gew_p1,2));
ef_gew_p1_lo=squeeze(prctile(ef_gew_p1,lo,2));
ef_gew_p1_hi=squeeze(prctile(ef_gew_p1,hi,2));
ef_gew_p2_m=squeeze(median(ef_gew_p2,2));
ef_gew_p2_lo=squeeze(prctile(ef_gew_p2,lo,2));
ef_gew_p2_hi=squeeze(prctile(ef_gew_p2,hi,2));



%% figure
% set(0,'DefaultPlotFontSize', 14);
set(0,'defaultAxesFontSize', 13);
set(0,'defaultTextFontSize', 14);

colth=[128 0 0]/255; % color of theoretical value
colm=[86 86 158]/255; % color of median value
colsh=[189 185 219]/255; %color of shades
figure(1); clf
ymax=[6 0.6 0.4];
ymin=[-0.7 -0.1 -0.25];
for it=1:length(jjv)
    %ymax(it)=1.02*(max(max([ef_gew_p1_hi(:,it) ef_gew_p2_hi(:,it) ef_hi(:,it)])));
    %ymin(it)=1.02*(min(min([ef_gew_p1_lo(:,it) ef_gew_p2_lo(:,it) ef_lo(:,it)])));
    
    subplot(3,4,(it-1)*4+1);
    h1=area(freq,[ef_gew_p1_lo(:,it) ef_gew_p1_m(:,it)-ef_gew_p1_lo(:,it) ef_gew_p1_hi(:,it)-ef_gew_p1_m(:,it)]); hold on
    set(h1(1),'FaceColor','w'); set(h1(2),'FaceColor',colsh); set(h1(3),'FaceColor',colsh);
    set(h1(1),'EdgeColor','k'); set(h1(2),'EdgeColor',colsh); set(h1(3),'EdgeColor',colsh); set(h1(3),'LineWidth',1.5);
    plot(freq,ef_gew_p1_m(:,it),'Color',colm,'LineWidth',1.5);
    plot(freq,f(:,it),'Color',colth,'LineWidth',1.5);
    xlim([0 fc/2]); ylim([ymin(it) ymax(it)]); 
    xlabel('Frequency (Hz)'); ylabel(['f_{' int2str(iiv(it)) '\rightarrow' int2str(jjv(it)) '}'],'FontSize',16);
    set(gca,'XTick',[0,10,40,50]);
    set(gca,'XTickLabel',{'0','10','40','50'});
    if it==1, title(['VAR, p=' int2str(p1)]); 
    hleg=legend('','est., dispersion','true','est., median','true','location','NorthWest');
    set(hleg,'FontSize',9); %set(hleg,'visible','off');
    legend boxoff;
    end
    
    subplot(3,4,(it-1)*4+2);
    h2=area(freq,[ef_gew_p2_lo(:,it) ef_gew_p2_m(:,it)-ef_gew_p2_lo(:,it) ef_gew_p2_hi(:,it)-ef_gew_p2_m(:,it)]); hold on
    set(h2(1),'FaceColor','w'); set(h2(2),'FaceColor',colsh); set(h2(3),'FaceColor',colsh);
    set(h2(1),'EdgeColor','k'); set(h2(2),'EdgeColor',colsh); set(h2(3),'EdgeColor','w'); set(h2(3),'LineWidth',1.5);
    plot(freq,ef_gew_p2_m(:,it),'Color',colm,'LineWidth',1.5);
    plot(freq,f(:,it),'Color',colth,'LineWidth',1.5);
    xlim([0 fc/2]); ylim([ymin(it) ymax(it)]); 
    xlabel('Frequency (Hz)'); ylabel(['f_{' int2str(iiv(it)) '\rightarrow' int2str(jjv(it)) '}'],'FontSize',16);
    if it==1, title(['VAR, p=' int2str(p2)]); end
    set(gca,'XTick',[0,10,40,50]);
    set(gca,'XTickLabel',{'0','10','40','50'});
    
    subplot(3,4,(it-1)*4+3);
    h3=area(freq,[ef_lo(:,it) ef_m(:,it)-ef_lo(:,it) ef_hi(:,it)-ef_m(:,it)]); hold on
    set(h3(1),'FaceColor','w'); set(h3(2),'FaceColor',colsh); set(h3(3),'FaceColor',colsh);
    set(h3(1),'EdgeColor','k'); set(h3(2),'EdgeColor',colsh); set(h3(3),'EdgeColor','w'); set(h3(3),'LineWidth',1.5);
    plot(freq,ef_m(:,it),'Color',colm,'LineWidth',1.5);
    plot(freq,f(:,it),'Color',colth,'LineWidth',1.5);
    xlim([0 fc/2]); ylim([ymin(it) ymax(it)]); 
    xlabel('Frequency (Hz)'); ylabel(['f_{' int2str(iiv(it)) '\rightarrow' int2str(jjv(it)) '}'],'FontSize',16);
    if it==1, title(['ISS, p=' int2str(p1)]); end
    set(gca,'XTick',[0,10,40,50]);
    set(gca,'XTickLabel',{'0','10','40','50'});
    
    subplot(3,4,(it-1)*4+4);%plot 1st realization
    hold on
    plot(freq,ef_gew_p1(:,1,it),'Color',[64 128 0]/255,'LineWidth',1.5);
    plot(freq,ef_gew_p2(:,1,it),'Color',[200 0 200]/255,'LineWidth',1.5);
    plot(freq,ef(:,1,it),'Color',[255 128 0]/255,'LineWidth',1.5);
    plot(freq,f(:,it),'Color',colth,'LineWidth',1.5);
    xlim([0 fc/2]); xlabel('Frequency (Hz)'); ylabel(['f_{' int2str(iiv(it)) '\rightarrow' int2str(jjv(it)) '}'],'FontSize',16);
    if it==1, title(['single run']); 
    hleg=legend('VAR, p=3','VAR, p=20','SS, p=3','true','location','NorthWest');
    set(hleg,'FontSize',9); %set(hleg,'visible','off');
    legend boxoff;
    end
    set(gca,'XTick',[0,10,40,50]);
    set(gca,'XTickLabel',{'0','10','40','50'});
end

%%%% use export_fig toolbox to export
% set(gcf, 'Position', [1 1 1400 900]);
% set(gcf, 'Color', 'w');
% saveas(gcf, 'test.png');
% export_fig test2.png -m2.5




