function [mic,mim]= roi_mim(Cohroi, subinds1, subinds2)

regu = 0.000001; 
[~,~,nfreq] = size(Cohroi);

ipcs = numel(subinds1);
jpcs = numel(subinds2);
        
for ifq = 1:nfreq
  cs_red=[];
  cs_red1 = Cohroi(subinds1,subinds1,ifq);
  cs_red2 = Cohroi(subinds1,subinds2,ifq);
  cs_red3 = Cohroi(subinds2,subinds2,ifq);

  caainv=inv(real(cs_red1)+regu*eye(ipcs)*mean(diag(real(cs_red1))));
  cab=imag(cs_red2);
  cbbinv=inv(real(cs_red3)+regu*eye(jpcs)*mean(diag(real(cs_red3))));
  X=cab*cbbinv*cab';
  % MIM Ewald Eq. 14
  mim(ifq)=(trace(caainv*X));
  caainvsqrt=sqrtm(caainv);
  Y=caainvsqrt*X*caainvsqrt; %Eq. 23
  [~,s,~]=svd(Y);
  % MIC
  mic(ifq)=sqrt(s(1,1));
end        

