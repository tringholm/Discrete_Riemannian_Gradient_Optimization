clear all
% close all
load('convtest17ApRadapt.mat')
Vhist2 = Vhist;
Vmin2 = Vmin;
load('convtest17ApR.mat')
Vhist = Vhist(1:10001);
n = 1:10:length(Vhist);
figure; loglog(n,(Vhist(1:10:end) - Vmin)./(Vhist(1) - Vmin),'b','linewidth',1.5)
hold on
loglog(1:length(Vhist2),(Vhist2 - Vmin2)./(Vhist2(1) - Vmin2),'b-.','linewidth',1.5)
loglog(n,1./n.^0.5,'k-.','linewidth',1.5)
loglog(n(1:70),1./n(1:70).^2,'k--','linewidth',1.5)

axis tight
legend('DRG','DRG-ADAPT','1/n^{1/2}','1/n^2','location','best')
xlabel('Iteration no.')
ylabel('Optimality error')
set(gcf,'paperposition',[0 0 5 4.2])
set(gca,'linewidth',1);
filename = ['InSARConv'];
print('-depsc','-r300',[filename '.eps'])