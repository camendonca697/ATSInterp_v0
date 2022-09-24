clear
load ftt.res -ascii
load tt.res -ascii
load gz2x.res -ascii
load x0.res -ascii
load zo.res -ascii
load pp.res -ascii
load zmax.res -ascii
%
[nx ny]=size(tt);
gz2=gz2x;
xi=min(x0);
xf=max(x0);
dd=mean(diff(x0));
%------------------------------------------------------
iFm0 = gz2 < 0.0;
raTm0=gz2(iFm0)./tt(iFm0);
hhTm0=sqrt(-1./raTm0)+zo;
ZapNaN=nan(nx,1);ZapNaN(iFm0)=hhTm0;
ZapZer=1000*ones(nx,1);ZapZer(iFm0)=hhTm0;
[pk,iloc,wpk,p]=findpeaks(-ZapZer);
%
figure
subplot(331)
plot(x0,tt,'-k',x0,ftt,'-r','LineWidth',1.5);
ww=axis;hold on;w(1:2)=[xi xf];plot(ww(1:2),[0 0],'--k');hold off
subplot(334)
plot(x0,gz2,'-k','LineWidth',1.5);
ww=[xi xf 5*floor(min(gz2)/5) 5*ceil(max(gz2)/5)];axis(ww)
hold on;
plot(ww(1:2),[0 0],'--k');hold off
subplot(337)
axis ij
hold on
ww=plotaprismaI(pp,'r',zmax);
plot(x0,ZapNaN,'-k','LineWidth',1.0)
plot(x0(iloc),ZapZer(iloc),'ok','MarkerFaceColor','w','MarkerSize',4)
hold off
axis image;axis([xi xf 0 zmax]);
%
subplot(334)
hold on
plot(x0(iloc),gz2(iloc),'ok','MarkerFaceColor','w','MarkerSize',4)
hold off
