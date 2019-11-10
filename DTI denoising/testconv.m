clear 
close all
addpath Eig3Folder
addpath example_dwi
addpath NIfTI_20140122
addpath fanDTasia
addpath export_fig

nii_file = load_untouch_nii('nifti_dt.nii');
data = nii_file.img;
clear nii_file;
k = 28; % 28 is standard case
slice = data(:,:,k,:,:);
slice2 = reshape(slice,112,112,6);

tensor = zeros(3,3,112,112);
mask = zeros(112);
for i = 1:112
    for j = 1:112
        tensor(1,1,i,j) = slice2(i,j,1);
        tensor(2,1,i,j) = slice2(i,j,2);
        tensor(1,2,i,j) = slice2(i,j,2);
        tensor(2,2,i,j) = slice2(i,j,3);
        tensor(3,1,i,j) = slice2(i,j,4);
        tensor(1,3,i,j) = slice2(i,j,4);
        tensor(3,2,i,j) = slice2(i,j,5);
        tensor(2,3,i,j) = slice2(i,j,5);
        tensor(3,3,i,j) = slice2(i,j,6);
        if any(any(tensor(:,:,i,j)))
            mask(i,j) = 1;
            [V,D] = eig(tensor(:,:,i,j));
            nmeans = 0;
            if any(any(D < 0))
%                 D(D < 0) = 1E-9;
%                 tensor(:,:,i,j) = V*D*V';
%                 tensor(:,:,i,j) = zeros(3);
                if i > 1
                    tensor(:,:,i,j) = tensor(:,:,i,j) + tensor(:,:,i-1,j);
                    nmeans = nmeans + 1;
                end
                if i < 112
                    tensor(:,:,i,j) = tensor(:,:,i,j) + tensor(:,:,i+1,j);
                    nmeans = nmeans + 1;
                end
                if j > 1
                    tensor(:,:,i,j) = tensor(:,:,i,j) + tensor(:,:,i,j-1);
                    nmeans = nmeans + 1;
                end
                if j < 112
                    tensor(:,:,i,j) = tensor(:,:,i,j) + tensor(:,:,i,j+1);
                    nmeans = nmeans + 1;
                end
                tensor(:,:,i,j) = tensor(:,:,i,j)/nmeans;
                [V,D] = eig(tensor(:,:,i,j));
                if any(any(D < 0))
                    D(D < 0) = 1E-9;
                end
                tensor(:,:,i,j) = V*D*V';
            end
        end
    end
end

yind = 34:73;
xind = 28:80;
yind = 45:50;
xind = 45:50;
% lowyind = find(sum(mask),1,'first');
% higyind = find(sum(mask),1,'last');
% lowxind = find(sum(mask,2),1,'first');
% higxind = find(sum(mask,2),1,'last');
% xind = lowxind:higxind;
% yind = lowyind:higyind;

tensor = tensor(:,:,xind,yind);
mask = mask(xind,yind);


tensor2 = tensor./max(max(max(max(tensor))));
tensor2 = tensor2*5;


gamma = 0.11;
% gamma = 0.05
dt = 0.001; % 0.01, 0.05, combi 
tol = 1E-5;
xtol = 1E-6;
p = 2;
q = 1;
g = tensor./max(max(max(max(tensor))));
[~,~,M,N] = size(g);
gg = zeros(6,M,N);
for i = 1:M
    for j = 1:N
       uT = g(:,:,i,j);
       if i == 48 && j == 22
           i
       end
       gg(:,i,j) = [uT(1); uT(2); uT(3); uT(5); uT(6); uT(9)];
    end 
end
g = gg;
% mex DGstep.c
% Vmin = 3.669165613289171e+02; % alpha = 0.11, full slice no. 28. Got with dt = 0.03
% Vmin = 2.125986539632652e+02
Vmin = 3.668990150487759e+02;
Vmin = 0;
Vminfact = 1E-8;
tic
[uu,Vhist2] = DTI_denoising10(g,mask,p,q,gamma,tol,xtol,dt);
% [uu, Vhist] = DTI_denoising_timing(g,mask,p,q,gamma,tol,xtol,dt,Vmin,Vminfact);
toc
Vhist = Vhist(Vhist>0);

% figure
% Vopt = 3.669185052822795e+02; % Probably more realistic, alpha = 0.11, full slice no. 28. Got with dt = 0.025
% Vopt = 3.668990150487759e+02; % alpha 0.11, full slice no. 28. Got with dt*k/k+1 and dt0 = 0.5
% loglog((Vhist-Vopt)./(Vhist(1) - Vopt),'linewidth',1.5)
% hold on
x = 1:(length(Vhist));
% loglog(x.^(-2),'linewidth',1.5)
% loglog(x.^(-1),'linewidth',1.5)
% axis tight
% legend('\tau = 0.05','1/n^2','1/n','\tau = 0.01','\tau = 0.05/0.01','location','northeast')
% xlabel('Iteration no.')
% ylabel('Optimality error')
% set(gcf,'paperposition',[0 0 5 4.2])
% set(gca,'linewidth',1);
% filename = ['DTIConv'];
% print('-depsc','-r300',[filename '.eps'])
% Vhist(100)-Vopt % 0.065 for dt = 0.025, rr = 1.1

u = zeros(3,3,M,N);
for i = 1:M
    for j = 1:N
       uT = uu(:,i,j);
       u(:,:,i,j) = [uT(1) uT(2) uT(3);
                     uT(2) uT(4) uT(5);
                     uT(3) uT(5) uT(6)];
    end 
end
Vhist(end)
length(Vhist)
% u = u./max(max(max(max(u))));
u = u*5;
figure
plotDTI(u)