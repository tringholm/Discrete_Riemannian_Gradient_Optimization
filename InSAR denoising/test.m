clear all
close all
[g, map] = imread('Vesuflat.gif');
g = map(g);
% g = g/max(max(g));
g = 2*pi*g - pi;
p = 2;
q = 1;
alpha = 0.3;
tol = 1E-8;
xtol = 1E-12;
dt = 0.25;
N = 40;
g = g(end-149:end,1:150);
% g = g(1:N,1:N);
plotnsave(g,'Original','InSAR_orig');
epsilon = 1E-2;

% tic
% [u,Vhist] = InSAR_DG(g,p,q,alpha,tol,xtol,dt,epsilon);
% toc
tic
% Vmin = 4.869743527989891e+03; %  Gotten with xtol = 1E-12
% Vmin = 4.866933345067895e+03; % alpha = 0.3, xtol = 1E-12;
% Vmin = 4.866723524364425e+03; % alpha = 0.3, xtol = 1E-12; dt = 0.001
% Vmin = 4.864920658967383e+03; % alpha = 0.3, xtol = 1E-12; dt = 0.002
% Vmin = 4.863957071621246e+03; % alpha = 0.3, xtol = 1E-12; dt = 0.002, 40 000 iterations
Vmin = 4.862372164846530e+03; % alpha = 0.3, xtol = 1E-12; dt = 0.05 (/2 every 200 iterations), 6637 iterations

% Vmin = 0
Vminfact = 1+1E-6;
dt = 0.05
[u,Vhist] = InSAR_DG_NoSmooth_timing(g,p,q,alpha,tol,xtol,dt,Vmin,Vminfact);
% [u,Vhist] = InSAR_DG_OLD(g,p,q,alpha,tol,dt);

toc
Vhist(end)
% addpath 'PhaseUnwrapping2D'
% g_uwr = phase_unwrap(g);
% hsv2 = imadjust(hsv,[0, 0, 0; 1, 1, 1],[0.25, 0.25, 0.25; 0.75, 0.75, 0.75]);
% plotnsave(g_uwr,'Unwrapped original','InSAR_orig_unw');
% figure; imagesc(u); title(['Reconstructed with L' int2str(p) ' fidelity term.']); colorbar;
% plotnsave(u,['Reconstructed with L' int2str(p) ' fidelity term.'],'InSAR_L2');
% u_uwr = phase_unwrap(u);
% plotnsave(u_uwr,['Unwrapped, reconstructed with L' int2str(p) ' fidelity term.'],'InSAR_L2_uwr');


figure; loglog((Vhist - Vmin)./(Vhist(1) - Vmin),'linewidth',1.5)
hold on; loglog(1./(1:length(Vhist)).^2,'linewidth',1.5)
hold on; loglog(1./(1:length(Vhist)),'linewidth',1.5)
axis tight
legend('DG','1/n^2','1/n','location','northeast')
xlabel('Iteration no.')
ylabel('Optimality error')
set(gcf,'paperposition',[0 0 5 4.2])
set(gca,'linewidth',1);
filename = ['InSARConv'];
print('-depsc','-r300',[filename '.eps'])