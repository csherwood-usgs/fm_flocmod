% fm_complot - Compare model results
load('fm.dat')
load('C:\\crs\\src\\flocmod\\flocmod.dat')
figure(1)
clf
subplot(211)
h1=plot(fm(:,1)./60,fm(:,2),'-k','linewidth',3)
set(h1,'color',[.5 .5 .5])
hold on
h2=plot(flocmod(:,1)./60,flocmod(:,2),'-k','linewidth',1)
ylabel('G (s^{-1})')
subplot(212)
h3=plot(fm(:,1)./60,fm(:,3),':g','linewidth',2)
set(h3,'color',[.7 .5 .5])
hold on
h4=plot(flocmod(:,1)./60,flocmod(:,3),'-r')
xlabel('Time (min)')
ylabel('D_{50} ({\mum})')
%%
subplot(212)
legend([h4;h3],'Fortran, double','Matlab')