% Denoise InSAR image 
 
function [u,Vhist] = InSAR_DG_OLD(g,p,q,gamma,tol,dt)
u = g;
% u = mod(g + randn(size(g))+pi,2*pi)-pi;
 
Tmax = 1000;
residual = inf;
 
Vhist = zeros(Tmax,1);
Vhist(1) = Veval(u,g,p,q,gamma);
tic
imfig = figure;
[M,N] = size(g);
uxm = 0; uym = 0; uxp = 0; uyp = 0;
 
dp = zeros(M,N);
dqx = zeros(M,N+1);
dqy = zeros(M+1,N);
for i = 1:M
   for j = 1:N
       dp(i,j) = d_ang(u(i,j),g(i,j))^p;
   end
end
for i = 1:M
   for j = 2:N
       dqx(i,j) = d_ang(u(i,j),u(i,j-1))^q;
   end
end
for i = 2:M
   for j = 1:N
       dqy(i,j) = d_ang(u(i,j),u(i-1,j))^q;
   end
end
 
for k = 1:Tmax
   if  residual < tol
      break
   end
   u_old = u; 
   for i = 1:M
       for j = 1:N
           uo = u(i,j);
           if i < M
           uyp = u(i+1,j);
           end
           if i > 1
           uym = u(i-1,j);
           end
           if j < N
           uxp = u(i,j+1);
           end
           if j > 1
           uxm = u(i,j-1);
           end
           Vo = 1/p*dp(i,j) + gamma/q*(dqx(i,j) + dqx(i,j+1) + dqy(i,j) + dqy(i+1,j));
           solveFxn = @(alp) alp + dt*Vdiff(alp,uo,uxp,uyp,uxm,uym,g(i,j),p,q,gamma,i,j,M,N,Vo)/alp;
           [alpha,fval] = fzero(solveFxn,0.01);%,optimset('Display','iter'));
%            alpha
           u(i,j) = mod(u(i,j) + alpha + pi,2*pi) - pi;
       end
   end
   if any(any(u > pi)) || any(any(u < -pi))
       disp('oooa')
   end
   Vhist(k+1) = Veval(u,g,p,q,gamma);
   Vhist(k+1)
   figure(imfig); imagesc(u); title('Intermediate plot');
   pause(0.01)
   residual = (Vhist(k) - Vhist(k+1))/Vhist(1);
end
toc
residual
Vhist = Vhist(~~Vhist);
end
 
function out = Vdiff(alp,u,uxp,uyp,uxm,uym,g,p,q,gamma,i,j,M,N,Vo)
un = mod(u + alp + pi,2*pi) - pi;
out = 1/p*d_ang(un,g)^p - 1/p*d_ang(u,g)^p;
if j < N
    out = out + gamma/q*(d_ang(un,uxp)^q - d_ang(u,uxp)^q);
end
if j > 1
    out = out + gamma/q*(d_ang(uxm,un)^q - d_ang(uxm,u)^q);
end
if i < M
    out = out + gamma/q*(d_ang(un,uyp)^q - d_ang(u,uyp)^q);
end
if i > 1
    out = out + gamma/q*(d_ang(uym,un)^q - d_ang(uym,u)^q);
end
% Vo - (1/p*d_ang(u,g)^p + gamma/q*(d_ang(u,uxp)^q + d_ang(uxm,u)^q + d_ang(u,uyp)^q + d_ang(uym,u)^q))
end
 
function out = d_ang(phi,theta)
out = pi - abs(sqrt((phi-theta)^2 + 1E-6) -pi); 
end
 
 
function out = Veval(u,g,p,q,gamma)
[M,N] = size(u);
dp = zeros(size(u));
dqx = zeros(size(u));
dqy = zeros(size(u));
for i = 1:M
   for j = 1:N
       dp(i,j) = d_ang(u(i,j),g(i,j))^p;
   end
end
for i = 1:(M-1)
   for j = 1:N
       dqx(i,j) = d_ang(u(i,j),u(i+1,j))^q;
   end
end
for i = 1:M
   for j = 1:(N-1)
       dqy(i,j) = d_ang(u(i,j),u(i,j+1))^q;
   end
end
out = 1/p*sum(sum(dp)) + gamma/q*(sum(sum(dqx)) + sum(sum(dqy)));
end