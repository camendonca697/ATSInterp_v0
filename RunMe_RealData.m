% code name
% ATSInterp version 0.0  - Automatic Thin-Sheet Interpretation
%
% manuscript
% THIN SHEET AUTOMATIC INTERPRETATION USING REGULARIZED SECOND DERIVATIVES   
% 
% REAL DATA APPLICATION CODE
%
clear
clc
help RunMe_RealData
% 1) input of profile data
load x0.dat -ascii
load ft.dat -ascii
load tt.dat -ascii
dd=50;   % data space along profile (m)
zo=-100; % sensor mean height with respect to the ground surface (m) 
nx=length(tt);
% data fitting from thin-sheet model in Cavalcante et al., 2020 paper
load JH_DataFitting.dat -ascii
%
% 2) evaluate regularized 2nd-derivative
aR1 = rdiff(tt(:),dd,'tikhonov',0.1);
gz2 = rdiff(aR1,dd,'tikhonov',0.01);
%
% 3) evaluate intervals with negative 2nd-derivative
iFm0 = gz2 < 0.0;
%
% 4)evaluate za function at intervals with negative 2nd-derivative
raTm0=tt(iFm0)./gz2(iFm0);
hhTm0=sqrt(-raTm0)+zo;
 A0m0=-[hhTm0.^3].*gz2(iFm0)/200;
%
% 5) let NaN dummies where the 2nd-derivative is non-negative
xNaN=nan(nx,1);
gz2NaN=xNaN;gz2NaN(iFm0)=gz2(iFm0);
ZapNaN=xNaN;ZapNaN(iFm0)=hhTm0;
 A0NaN=xNaN; A0NaN(iFm0)=A0m0;
%
% 6) evaluate local minima in the 2nd-derivative profile
[pk,iloc,wpk,p]=findpeaks(-gz2NaN);
pZa(:,1:3)=[x0(iloc) ZapNaN(iloc) A0NaN(iloc)];
%
% 7) xls output with (x0,z0,a0) estimates
xlswrite('ZaResults.xls',pZa);
%
% 8) composed figure as in the manuscript
xi=min(x0)/1000;
xf=max(x0)/1000;
figure
subplot(311)
plot(x0(1:4:end)/1000,ft(1:4:end),'sk','MarkerSize',4);
hold on
plot(JH_DataFitting(:,1)/1000,JH_DataFitting(:,2),'-r','LineWidth',1.5);hold off
legend('observed','evaluated');
axis([xi xf -300 600])
ylabel('TFA (nT)')
subplot(613)
plot(x0/1000,tt,'-k','LineWidth',1.5);
hold on;plot(x0(iloc)/1000,tt(iloc),'ok','MarkerFaceColor','w','MarkerSize',4);hold off
axis([xi xf 0 max(tt)])
ylabel('AMA (nT)')
subplot(614)
plot(x0/1000,gz2,'-k','LineWidth',1.5);
ww=[xi xf min(gz2) max(gz2)];hold on;
plot(x0(iloc)/1000,gz2(iloc),'ok','MarkerFaceColor','w','MarkerSize',4)
plot([xi xf],[0 0],'--k');hold off  
axis([xi xf min(gz2) max(gz2)])
ylabel('(nT/km^2)') 
subplot(313)
axis ij
hold on
zmax=1000;
aux01_ModelJH;% model as in Cavalcante et al.,2020 paper
plot(x0/1000,ZapNaN,'-k','LineWidth',1.5)
plot(x0(iloc)/1000,ZapNaN(iloc),'ok','MarkerFaceColor','w','MarkerSize',4)
for k=10:10:50;text(ps(k,3)/1000-0.3,ps(k,4)-0.3,num2str(k));end
axis([xi xf 0 zmax]);
hold off
pause(4)
%
% 9) close up look to local-minimum points
figure
subplot(211)
plot(x0/1000,tt,'-k','LineWidth',1.5);hold on
plot(x0(iloc)/1000,tt(iloc),'ok','MarkerFaceColor','y','MarkerSize',4);hold off
axis([xi xf 0 max(tt)]);grid on
text(3,500,['Local peak points (yellow): ' num2str(length(iloc))])
ylabel('AMA (nT)')
subplot(212)
plot(x0/1000,gz2,'-k','LineWidth',1.5);
ww=[xi xf min(gz2) max(gz2)];hold on;
plot(x0(iloc)/1000,gz2(iloc),'ok','MarkerFaceColor','y','MarkerSize',4)
plot([xi xf],[0 0],'--k');hold off  
axis([xi xf min(gz2) max(gz2)]);grid on
text(17,0.0026,['Local minimum points (yellow): ' num2str(length(iloc))])
ylabel('2nd derivative (nT/km^2)') 
