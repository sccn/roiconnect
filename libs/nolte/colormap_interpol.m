function newmap=colormap_interpol(map,fak);
% creates a finer colormap consisting of fak times more colors)

[n,ndum]=size(map);

newmap(1,:)=map(1,:);

m=1;
for i=2:n;
    for k=1:fak
        m=m+1;
        newmap(m,:)=((fak-k)*map(i-1,:)+k*map(i,:))/fak;
    end
end

return;