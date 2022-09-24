clear
% X Y T |T|
load Perfil_JH.XYZ
ax=Perfil_JH(:,1);
ay=Perfil_JH(:,2);
ft=Perfil_JH(:,3);
tt=Perfil_JH(:,4);
ax=ax-ax(1,1);
ay=ay-ay(1,1);
nx=length(ax(:,1));
x0=zeros(nx,1);
for i=2:nx
    x0(i,1)=x0(i-1)+sqrt([ax(i)-ax(i-1)]^2+[ay(i)-ay(i-1)]^2);
end
xi=min(x0)/1000;
xf=max(x0)/1000;
% calculado # POS (m) T (nT) |T|(nT) i_m (graus)
load calculado_JH.dat -ascii
x0c=calculado_JH(:,1);
ftc=calculado_JH(:,2);
ttc=calculado_JH(:,3);
%
subplot(211)
plot(x0/1000,ft,'sk','MarkerSize',2);hold on
plot(calculado_JH(:,1)/1000,calculado_JH(:,2),'-r');hold off
legend('obs','calculado');
axis([xi xf -500 500])
xlabel(' x (km)');
ylabel('T_t (nT)')
%
% MODELO COM OS DIQUES (ARQUIVO EXITMOD_JH)
% a (A) i_m (graus) x0 (m) z0 (m)
load exitmod_JH.dat -ascii
ps=exitmod_JH;
ns=length(ps);
subplot(212)
zmax=1000;
plot([ps(1,3) ps(1,3)]/1000,[ps(1,4) zmax],'-k');
hold on
for i=2:ns
    ch='b-';if ps(i,2)>0;ch='r-';end
    plot([ps(i,3) ps(i,3)]/1000,[ps(i,4) zmax],ch,'LineWidth',1.5);
end
axis ij
xlabel(' x (km)')
ylabel('z (m)')
axis([xi xf 0 zmax])
%
