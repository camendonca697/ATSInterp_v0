% para B=30000 nT e KT=60 temos:
%                 X=6x10(-3)SI
%                 Mi= 0.1432 A/m  
%                 60  --> 0.1432
%                  kt  --> Mi
%                  kt = (60*Mi)/0.1432 = 418.8790 Mi 
clear
xi=-40;xf=40;dd=0.25;
x0=[xi:dd:xf]';
nx=length(x0);
z0=zeros(nx,1);
HH=230;
zT=10*dd;
b=3*dd;
zmax=10;
%
Tinc=30;
MincI=Tinc;
MincR=90;
zo=0;
%
clc
pp=[0  zT  b   70    1   MincI  0 Tinc HH];
AV0=pp(3)*2*pp(5);
pp(:,5)=418.8790*pp(:,5); % conversao A/m to kt
[H Z ftt]=multiprism(x0,pp);
TT=sqrt(H.*H+Z.*Z);
%
aR1 = rdiff(TT,dd,'tikhonov',0.001);
gz2 = rdiff(aR1,dd,'tikhonov',0.001);
iFm0 = gz2 < 0.0;
raTm0=TT(iFm0)./gz2(iFm0);
hhTm0=sqrt(-raTm0)+zo;
X0=NaN(nx,1);Y0=X0;X0(iFm0,1)=x0(iFm0);Y0(iFm0,1)=hhTm0;
[aa ia]=min(abs(x0));
[zh izh]=min(Y0);
A0=-[zh^2]*gz2(ia)/100;
pT=[A0 AV0];
disp(['A0=' num2str(pT(1,1)) '  V0=' num2str(pT(1,2))]);
XI=-40*dd;XF=-XI;
%
xxi=-zh/sqrt(2);xxf=-xxi;xx=[xxi+dd/2:dd/2:xxf-dd/2]';xx2=xx.*xx;
a1=(xx.^2+zh*zh);
a2=sqrt(zh*zh-2*xx2);
zap=a1./a2;
%
figure
subplot(621);plot(x0,TT ,'-k');ww=axis;ww(1:2)=[XI XF];axis(ww)
subplot(623);plot(x0,gz2,'-k',[XI XF],[0 0],'--k');hold on
plot(x0(ia),gz2(ia),'ok','MarkerSize',6,'MarkerFaceColor','b');hold off
ww=axis;ww(1:2)=[XI XF];axis(ww)
subplot(325);
axis ij
hold on
ww=plotaprismaI(pp,'y',zmax);
plot(xx,zap,'-c','LineWidth',4)
plot(X0,Y0,'-k','LineWidth',2.0)
legend(' ',' ',' ')
plot([0 0],[zh zmax],'-b','LineWidth',2)
plot(X0(izh),Y0(izh),'ok','MarkerSize',6,'MarkerFaceColor','b')
hold off
axis image;axis([XI XF 0 zmax]);

