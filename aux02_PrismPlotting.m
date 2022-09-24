function v = PrismPlotting(pp,ch,zmax)
%pp=[1x0 2h 3b 4delta 5KT 6Minc 7alfa]
% para B=30000 nT e KT=60 temos:
%                 X=6x10(-3)SI
%                 Mi= 0.1432 A/m  
%                 60  --> 0.1432
%                 KT  --> Mi
%                 KT = (60*Mi)/0.1432 = 418.8790 Mi 
%                 Mi = KT/418.8790  
np=length(pp(:,1));
for k=1:np
    zm=min(zmax,pp(k,9));
    x0=pp(k,1);
    h=pp(k,2);
    b=pp(k,3);
   KT=pp(k,5);
    v(k,1:3)=[x0 h 4*pi*2*b*KT/418.8790]; 
    delta=pp(k,4);
    wx=[x0-b x0+b x0+b+cosd(delta)*zm x0-b+cosd(delta)*zm]';
    wz=[h h zm zm]';
    fill(wx,wz,ch)
end
axis ij
end

