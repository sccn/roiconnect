%% stima le funzioni f e fcond col formalismo di Geweke (1982, 1984)

% Y: Q*N data
% Mv: number of series in each block
% p: order of the MVAR model to be used to fit the data
% ii, jj: indexes of the blocks to which estimate the measure: from Yii to Yjj

function [gewF,gewFc,f] = Geweke_f(Y,Mv,p,ii,jj,nfft,fc)

Q=size(Y,1);
M=length(Mv);


%% estimated f Geweke bivariate measure
i1=sum(Mv(1:ii)); i0=i1-Mv(ii)+1; %indici dei due blocchi di interesse
j1=sum(Mv(1:jj)); j0=j1-Mv(jj)+1;
Yb=Y([i0:i1 j0:j1],:); % solo i segnali dei 2 blocchi di interesse
Mvb=[Mv(ii) Mv(jj)]'; % vettore bivariato con le dimensioni dei 2 blocchi di interesse
[Amtmp,Sutmp]=idMVAR(Yb,p,0); %stimo modello MVAR solo sui 2 blocchi di interesse
[bDCtmp,bPDCtmp,emF2,mGtmp,bStmp,bPtmp,bHtmp,bAftmp,f] = block_fdMVAR(Amtmp,Sutmp,Mvb,nfft,fc); %tiro fuori emF2
gewF=squeeze(emF2(2,1,:)); %% misura bivariata da Yii a Yjj !!!

%% estimated f Geweke CONDITIONAL bivariate measure (è trivariate in realtà)
%%% NB: questa calcola da Yi a Yj !!!
ki=[]; %altri indici (non Yii, non Yjj)
for k=1:Q
    if sum(k==[i0:i1 j0:j1])==0, ki=[ki k]; end
end
Yb2=Y([j0:j1 ki],:); %tutti i segnali tranne quelli di Yi
Mvb2=[Mv(jj) Q-Mv(jj)-Mv(ii)]'; 
[Amtmp,Sutmp]=idMVAR(Yb2,p,0);
[bDCtmp,bPDCtmp,mFtmp,mGtmp,bStmp,bPtmp,GG] = block_fdMVAR(Amtmp,Sutmp,Mvb2,nfft,fc); % primo giro: G
sigmateta=Sutmp(1:Mv(jj),1:Mv(jj));

Yc=Y([j0:j1 i0:i1 ki],:);
Mvc=[Mv(jj) Mv(ii) Q-Mv(jj)-Mv(ii)]';
[Amtmp,Sutmp]=idMVAR(Yc,p,0);
[bDCtmp,bPDCtmp,mFtmp,mGtmp,bStmp,bPtmp,HH] = block_fdMVAR(Amtmp,Sutmp,Mvc,nfft,fc); % secondo giro: H
sigmaUx=Sutmp(1:Mv(jj),1:Mv(jj));

gewFc=NaN*ones(nfft,1); % Geweke conditional F da Yii a Yjj
for n=1:nfft
    Gtot=[GG{1,1,n} zeros(Mvc(1),Mvc(2)) GG{1,2,n}; zeros(Mvc(2),Mvc(1)) eye(Mvc(2)) zeros(Mvc(2),Mvc(3)); GG{2,1,n} zeros(Mvc(3),Mvc(2)) GG{2,2,n}];
    Htot=[HH{1,1,n} HH{1,2,n} HH{1,3,n}; HH{2,1,n} HH{2,2,n} HH{2,3,n}; HH{3,1,n} HH{3,2,n} HH{3,3,n}];
    QQ= inv(Gtot)*Htot;
    gewFc(n)=log( det(sigmateta) / abs( det(QQ(1:Mvc(1),1:Mvc(1))*sigmaUx* QQ(1:Mvc(1),1:Mvc(1))') ) );
end

