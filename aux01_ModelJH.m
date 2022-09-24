% thin sheet model as in Figure 10a of
% Cavalcante et al., Journal of Hydrology, 2020, DOI:10.1016/j.jhydrol.2020.125079.
%
% a (A) i_m (graus) x0 (m) z0 (m)
load JH_model.dat -ascii
ps=JH_model;
ns=length(ps);
plot([ps(1,3) ps(1,3)]/1000,[ps(1,4) zmax],'-k');
hold on
for i=2:ns
    ch='-r';if ps(i,2)>0;ch='b-';end
    plot([ps(i,3) ps(i,3)]/1000,[ps(i,4) zmax],ch,'LineWidth',2);
end
axis ij
xlabel('Distance (km)')
ylabel('Depth (m)')
axis([xi xf 0 zmax])
