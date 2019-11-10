function [u,Vhist] = DTI_denoising_mex3(g,mask,p,q,gamma,tol,xtol,dt)
u = zeros(size(g));
u_old = zeros(size(g));
u_old(:) = g(:);
u(:) = g(:);

[~,M,N] = size(g);

dp = zeros(M,N);
dqx = zeros(M,N+1);
dqy = zeros(M+1,N);
Suinvprod = zeros(6,6,M,N);
Sginvprod = zeros(6,6,M,N);


for i = 1:M
    for j = 1:N
        if mask(i,j)
            uT = u(:,i,j);
            uT = [uT(1) uT(2) uT(3);
                  uT(2) uT(4) uT(5);
                  uT(3) uT(5) uT(6)];
            Suinvprod(:,:,i,j) = prodform(inv(sqrtm(uT)));
            Sginvprod(:,:,i,j) = Suinvprod(:,:,i,j);
        end
    end
end
for i = 1:M
    for j = 1:N
        if mask(i,j)
            if ~any(any(isreal(Suinvprod(:,:,i,j))))
                i
                j
            end
            dp(i,j) = 1/p*d_M(Suinvprod(:,:,i,j),g(:,i,j))^p;
        end
    end
end
for i = 1:M
    for j = 2:N
        if mask(i,j) && mask(i,j-1)
            dqx(i,j) = d_M(Suinvprod(:,:,i,j),u(:,i,j-1))^q;
        end
    end
end
for i = 2:M
    for j = 1:N
        if mask(i,j) && mask(i-1,j)
            dqy(i,j) = d_M(Suinvprod(:,:,i,j),u(:,i-1,j))^q;
        end
    end
end
Tmax = 10000;

Vhist = zeros(Tmax,1);
Vhist(1) = 1/p*sum(sum(dp)) + gamma*(sum(sum(dqx)) + sum(sum(dqy)));
residual = inf;
prevres = residual;
rr = 1.1;
for k = 1:Tmax
%     tic
    if  residual < tol
        break
    end
    if k > 12
        dt = 0.01;
    end
    
    u = DGstep2(g,u,mask,dp,dqx,dqy,Suinvprod,Sginvprod,p,gamma,xtol,dt,M,N);
    Vhist(k+1) = sum(sum(dp)) + gamma/q*(sum(sum(dqx)) + sum(sum(dqy)));
    residual = (Vhist(k) - Vhist(k+1))/Vhist(1)
%     if prevres > residual
%         dt = dt*rr;        
%     else
%         dt = dt/rr^2;
%     end
%     dt

%     dt = dt*k/(k+1);
%     mean(mean(mean(abs(u - u_old))))
%     median(median(median(abs(u - u_old))))
%     max(max(max(abs(u - u_old))))
    u_old(:) = u(:);
    prevres = residual;
%     Vhist(k+1)
%     toc
end
residual
Vhist = Vhist(~~Vhist);
end



function out = d_M(XX,A)
A1 = XX(1,1)*A(1) + XX(2,1)*A(2) + XX(3,1)*A(3) + XX(4,1)*A(4) + XX(5,1)*A(5) + XX(6,1)*A(6);
A2 = XX(1,2)*A(1) + XX(2,2)*A(2) + XX(3,2)*A(3) + XX(4,2)*A(4) + XX(5,2)*A(5) + XX(6,2)*A(6);
A3 = XX(1,3)*A(1) + XX(2,3)*A(2) + XX(3,3)*A(3) + XX(4,3)*A(4) + XX(5,3)*A(5) + XX(6,3)*A(6);
A4 = XX(1,4)*A(1) + XX(2,4)*A(2) + XX(3,4)*A(3) + XX(4,4)*A(4) + XX(5,4)*A(5) + XX(6,4)*A(6);
A5 = XX(1,5)*A(1) + XX(2,5)*A(2) + XX(3,5)*A(3) + XX(4,5)*A(4) + XX(5,5)*A(5) + XX(6,5)*A(6);
A6 = XX(1,6)*A(1) + XX(2,6)*A(2) + XX(3,6)*A(3) + XX(4,6)*A(4) + XX(5,6)*A(5) + XX(6,6)*A(6);

p1 = A2^2 + A3^2 + A5^2;
q = (A1 + A4 + A6)/3;

A1 = A1 - q;
A4 = A4 - q;
A6 = A6 - q;

p2 = A1^2 + A4^2 + A6^2 + 2*p1;
p = 2*realsqrt(p2 / 6);

% Determinant
d = A5*(A1*A5 - 2*A2*A3) + A6*(A2^2 - A1*A4) + A4*A3^2;
r = -4*d/(p*p*p);

% In exact arithmetic for a symmetric matrix  -1 <= r <= 1
% but computation error can leave it slightly outside this range.
if (r < -1)
    eig1 = q + 0.5*p;
    eig3 = q - p;
elseif (r > 1)
    eig1 = q + p;
    eig3 = q - 0.5*p;
else
    phi = acos(r)/3;
    eig1 = q + p*cos(phi);
    eig3 = q + p*cos(phi + 2*pi/3);
end

eig2 = 3*q - eig1 - eig3;

out = sqrt(log(eig1)^2 + log(eig2)^2 + log(eig3)^2);

end

function XX = prodform(X)
x1 = X(1);
x2 = X(2);
x3 = X(3);
x4 = X(5);
x5 = X(6);
x6 = X(9);

x1x2 = x1*x2;
x1x3 = x1*x3;

x2x2 = x2*x2;
x2x3 = x2*x3;
x2x4 = x2*x4;
x2x5 = x2*x5;

x3x3 = x3*x3;
x3x5 = x3*x5;
x3x6 = x3*x6;

x4x5 = x4*x5;

x5x5 = x5*x5;
x5x6 = x5*x6;

x1x5 = x1*x5 + x2x3;
x3x4 = x3*x4 + x2x5;
x2x6 = x2*x6 + x3x5;
 
XX = [x1*x1,   x1x2,         x1x3,       x2x2,   x2x3,         x3x3;
      2*x1x2,  x1*x4 + x2x2, x1x5,       2*x2x4, x3x4,         2*x3x5;
      2*x1x3,  x1x5,         x1*x6+x3x3, 2*x2x5, x2x6,         2*x3x6;
      x2x2,    x2x4,         x2x5,       x4*x4,  x4x5,         x5x5;
      2*x2x3,  x3x4,         x2x6,       2*x4x5, x4*x6 + x5x5, 2*x5x6;
      x3x3,    x3x5,         x3x6,       x5x5,   x5x6,         x6*x6];
end