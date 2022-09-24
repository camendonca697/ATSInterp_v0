%
% Figure making for processing results
%
% 1) za-function for NOISE-FREE DATA
CF=200; % term 2pi/mu0 in equation 7 for fields in nT 
iFm0 = gz2 < 0.0; % intervals with negative 2nd-derivative
raTm0=TT(iFm0)./gz2(iFm0);
hhTm0=sqrt(-raTm0)+zo; % equation 5
A0m0=-[hhTm0.^3].*gz2(iFm0)/CF; % equation 7
xNaN=nan(nx,1); % dummies for plotting
gz2NaN=xNaN;gz2NaN(iFm0)= gz2(iFm0);
ZapNaN=xNaN;ZapNaN(iFm0)= hhTm0;
A0NaN =xNaN; A0NaN(iFm0)= A0m0;
[pk,iloc,wpk,p]=findpeaks(-gz2NaN);
pZa=[x0(iloc) ZapNaN(iloc) A0NaN(iloc)];
if nMOD==6;xlswrite('ZaResults.xls',pZa);end
%
% 2) za-function for NOISY DATA
iFm0R = gz2R < 0.0;
raTm0R=TTr(iFm0R)./gz2R(iFm0R);
hhTm0R=sqrt(-raTm0R)+zo;
A0m0R=-[hhTm0R.^3].*gz2(iFm0R)/CF;
ZapNaNR=xNaN;ZapNaNR(iFm0R)= hhTm0R;
A0NaNR=xNaN; A0NaNR(iFm0R) = A0m0R;
[pkR,ilocR,wpkR,pR]=findpeaks(-ZapNaNR);
%------------------------------------------------------
figure
subplot(421)
plot(x0,TTr,'-sc','LineWidth',1.0,'MarkerSize',3);
hold on
plot(x0,TT,'-k',x0,ftt,'-r','LineWidth',1.5)
hold off
legend('Noisy-AMA','AMA','TFA')
ylabel('Anomaly (nT)')
title(chm)
subplot(423)
plot(x0,gz2,'-k','LineWidth',1.5);
hold on;
plot(x0,gz2R,'-b','LineWidth',1.5);legend('true','regul.')
plot(x0(iloc),gz2(iloc),'ok','MarkerFaceColor','w','MarkerSize',4)
plot(x0(ilocR),gz2R(ilocR),'ob','MarkerFaceColor','w','MarkerSize',4)
hold off
ylabel('2nd drv (nT/km^2)')
subplot(425)
axis ij
hold on
pR=aux02_PrismPlotting(pp,'y',zmax);
plot(x0,ZapNaN,'-k','LineWidth',1.5);
axis([xi xf 0 zmax]);
plot(x0(iloc),ZapNaN(iloc),'ok','MarkerFaceColor','w','MarkerSize',4)
text(78.5,1.6,'true Z_a ')
hold off
ylabel('Depth (km)')
subplot(427)
axis ij
hold on
pR=aux02_PrismPlotting(pp,'y',zmax);
plot(x0,ZapNaNR,'-b','LineWidth',1.5)
axis([xi xf 0 zmax]);
plot(x0(ilocR),ZapNaN(ilocR),'ob','MarkerFaceColor','w','MarkerSize',4)
text(78.5,1.6,'regularized Z_a ')
hold off
ylabel('Depth (km)');xlabel('Distance (km)')
%
subplot(222)
plot(x0,gz2N,'-c',x0,gz2,'-k',x0,gz2R,'-b','LineWidth',1.0);
legend('unregularized','true','regularized')
ylabel('AMA 2nd derv. (nT/km2)')