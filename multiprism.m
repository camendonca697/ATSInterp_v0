function [H Z Tt]=multi_prism(x,PP)
% pp parameters for 2D prisms
% pp=[1x0 2h 3b 4delta 5KT 6Minc 7alfa 8Tinc 9HH] 
% x0 prism position along the profile
%  h depth to the top
%  b half width
% delta dip angle (degree)
% KT susceptibility (adm in SI) to geomagnetic field (nT) product
% Minc inclination of the magnetization
% alfa strike direction (clockwise counting from the North)
% Tinc inclination of the local main field
%  [distances in km, angles in degree]
% output
% H horizontal component of the magnetic field
% Z vertical component
% Tt total field anomaly (along direction Tinc)
np=length(PP(:,1));
nx=length(x(:,1));
H=zeros(nx,1);Z=H;
for k=1:np
    pp=PP(k,:);
    [fht fzt]=one_prism(x,nx,pp);
    pp(1)=pp(1)+pp(9)*cosd(pp(4));
    pp(2)=pp(2)+pp(9);
    [fhb fzb]=one_prism(x,nx,pp);
    H=H+fht-fhb;
    Z=Z+fzt-fzb;
end
Tt=cosd(pp(8))*H+sind(pp(8))*Z;
return
%
function  [H Z] = one_prism(x,nx,pp)
H=zeros(nx,1);Z=H;
wx=x-pp(1);
h=pp(2);
b=pp(3);
I=atand(tand(pp(6))/cosd(pp(7)));
teta=pp(6)-pp(4);
an1=cosd(pp(6));
an2=sind(pp(7));
C=2*pp(5)*sind(pp(4))*sqrt(1-an1*an1*an2*an2);
cteta=cosd(teta);
steta=sind(teta);
calfa=cosd(pp(7));
for i=1:nx
    a1=(wx(i)+b)/h;
    a2=(wx(i)-b)/h;
    b1=(wx(i)-b)^2+h*h;
    b2=(wx(i)+b)^2+h*h;
    wtg=[atan(a1)-atan(a2)];
    wln=log(b1/b2);
    Z(i,1)=C*(cteta*wtg-0.5*steta*wln);
    H(i,1)=C*(steta*wtg+0.5*cteta*wln);
end
return

