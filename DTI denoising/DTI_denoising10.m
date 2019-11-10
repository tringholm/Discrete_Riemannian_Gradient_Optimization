function [uu,Vhist] = DTI_denoising10(g,mask,p,q,gamma,tol,xtol,dt)
uu = g;
Tmax = 1000;
residual = inf;

uxm = 0; uym = 0; uxp = 0; uyp = 0;

[~,M,N] = size(g);

dp = zeros(M,N);
dqx = zeros(M,N+1);
dqy = zeros(M+1,N);
Suinvprod = zeros(6,6,M,N);
Sginvprod = zeros(6,6,M,N);


for i = 1:M
    for j = 1:N
        if mask(i,j)
            uT = uu(:,i,j);
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
            dp(i,j) = 1/p*d_M(Suinvprod(:,:,i,j),g(:,i,j))^p;
        end
    end
end
for i = 1:M
    for j = 2:N
        if mask(i,j) && mask(i,j-1)
            dqx(i,j) = d_M(Suinvprod(:,:,i,j),uu(:,i,j-1))^q;
        end
    end
end
for i = 2:M
    for j = 1:N
        if mask(i,j) && mask(i-1,j)
            dqy(i,j) = d_M(Suinvprod(:,:,i,j),uu(:,i-1,j))^q;
        end
    end
end
t = 0;
Vhist = zeros(Tmax,1);
Vhist(1) = 1/p*sum(sum(dp)) + gamma/q*(sum(sum(dqx)) + sum(sum(dqy)));
VV = Vhist(1);

% imfig = figure;
params = [p, gamma, 0, 0, dt];
masks = [0 0 0 0];

for k = 1:Tmax
    if  residual < tol
%         break
    end
%     tic
    for i = 1:M
%         i
        for j = 1:N
            if mask(i,j)
                masks(1) = 0; masks(2) = 0; masks(3) = 0; masks(4) = 0;
                uT = uu(:,i,j);
                gg = Sginvprod(:,:,i,j);
                if i < M
                    if mask(i+1,j)
                        uyp = Suinvprod(:,:,i+1,j);
                        masks(1) = 1;
                    end
                end
                if i > 1
                    if mask(i-1,j)
                        uym = Suinvprod(:,:,i-1,j);
                        masks(2) = 1;
                    end
                end
                if j < N
                    if mask(i,j+1) 
                        uxp = Suinvprod(:,:,i,j+1);
                        masks(3) = 1;
                    end
                end
                if j > 1
                    if mask(i,j-1)
                        uxm = Suinvprod(:,:,i,j-1);
                        masks(4) = 1;
                    end
                end
                if j == 5
                    pause(0.01);
                end
                for l = 1:6
                    params(3) = l;
                    uTinv = inv3x3(uT);
                    params(4) = dp(i,j) + gamma/q*(dqx(i,j) + dqx(i,j+1) + dqy(i,j) + dqy(i+1,j));
                    [alpha,  uT, dug, dxp, dxm, dyp, dym] = fzero_noExtra_spec([-dt dt],xtol,uT,uTinv,uxp,uyp,uxm,uym,gg,params,masks);
%                     uT
                    dp(i,j) = dug;
%                     VV = VV - alpha^2/dt;
                    if j > 1 && masks(4)
                        dqx(i,j) = dxm;
                    end
                    if j < N && masks(3)
                        dqx(i,j+1) = dxp;
                    end
                    if i > 1 && masks(2)
                        dqy(i,j) =  dym;
                    end
                    if i < M && masks(1)
                        dqy(i+1,j) =  dyp;
                    end
                end
                [V, L] = eignvec3x3S(uT);
                
                uu(:,i,j) = uT;
                Suinvprod(:,:,i,j) = prodform(V*(L.*V).');
                
                %             t = t + toc;
            end
        end
    end
%     Vhist(k+1) = VV;
    Vhist(k+1) = sum(sum(dp)) + gamma/q*(sum(sum(dqx)) + sum(sum(dqy)));
%     Vhist(k+1)
%     dt = dt*k/(k+1);
%     Vhist(k+1) - VV
%     uu = u*5;
%     clf(imfig);figure(imfig); plotDTI(uu); title('Intermediate plot');
%     drawnow
    residual = (Vhist(k) - Vhist(k+1))/Vhist(1);
    if ~mod(k,40)
        Vhist(k)
    end
%     toc
end
% t
residual
Vhist = Vhist(~~Vhist);
end

function uT = retraction(uT,uTinv,alpha,l)
C = 0.5*alpha^2;
if l == 1
        uT(1) = uT(1) + alpha + C*uTinv(1);
elseif l == 2
        uT(4) = uT(4) + alpha + C*uTinv(4);
elseif l == 3
        uT(6) = uT(6) + alpha + C*uTinv(6);
elseif l == 4
        uT(1) = uT(1) + C*uTinv(4);
        uT(2) = uT(2) + alpha + C*uTinv(2);
        uT(4) = uT(4) + C*uTinv(1);
elseif l == 5
        uT(4) = uT(4) + C*uTinv(6);
        uT(5) = uT(5) + alpha + C*uTinv(5); 
        uT(6) = uT(6) + C*uTinv(4);
else
        uT(1) = uT(1) + C*uTinv(6);
        uT(3) = uT(3) + alpha + C*uTinv(3); 
        uT(6) = uT(6) + C*uTinv(1);
end
end

% params: [p gamma l Vo dt]
function [out, uT, dug, dxp, dxm, dyp, dym] = Vdiff_solv(alpha,uT,uTinv,uxp,uyp,uxm,uym,g,params,masks)
l = params(3);
uT = retraction(uT,uTinv,alpha,l);

dug = 1/params(1)*d_M(g,uT)^params(1);
dyp = 0;
dym = 0;
dxp = 0;
dxm = 0;
if masks(1)
    dyp = d_M(uyp,uT);
end
if masks(2)
    dym = d_M(uym,uT);
end
if masks(3)
    dxp = d_M(uxp,uT);
end
if masks(4)
    dxm = d_M(uxm,uT);
end

out = dug - params(4) + params(2)*(dxp + dxm + dyp + dym);
out = alpha + params(5)*out/alpha;
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

function A = inv3x3(A)
a1 = A(1);
a2 = A(2);
a3 = A(3);
a4 = A(4);
a5 = A(5);
a6 = A(6);
aa1 = a6*a4-a5*a5;
aa2 = a5*a3-a2*a6;
aa3 = a5*a2-a4*a3;

d = a1*aa1 + a2*aa2 + a3*aa3;

A(1) = aa1;
A(2) = aa2;
A(3) = aa3;
A(4) = a6*a1 - a3*a3;
A(5) = a3*a2 - a5*a1;
A(6) = a4*a1 - a2*a2;
A = 1/d*A;
end

function [b, un, dug, dxp, dxm, dyp, dym] = fzero_noExtra_spec(x,tol,uT,uTinv,uxp,uyp,uxm,uym,g,params,masks)
badint = 0;
% Interval input
if (numel(x) == 2) 
 
    a = x(1);
    b = x(2);    
    [fa, una, duga, dxpa, dxma, dypa, dyma] = Vdiff_solv(a,uT,uTinv,uxp,uyp,uxm,uym,g,params,masks);
  
    [fb, un, dug, dxp, dxm, dyp, dym] = Vdiff_solv(b,uT,uTinv,uxp,uyp,uxm,uym,g,params,masks);
    
    if fa*fb > 0
%         disp('Wrong start interval');
        badint = 1;
        x = 0.01;
    end
    
    if ~fa
        b = a;
        un = una; dug = duga; dxp = dxpa; dxm = dxma; dyp = dypa; dym = dyma;
        return
    elseif ~fb
        % b = b;

        return
    end
end
    % Starting guess scalar input
if (numel(x) == 1) || badint

    
    [fx, un, dug, dxp, dxm, dyp, dym] = Vdiff_solv(x,uT,uTinv,uxp,uyp,uxm,uym,g,params,masks);

    if fx == 0
        b = x;
 
        return

    end
    
    if x ~= 0
        dx = x/50;
    else 
        dx = 1/50;
    end
    
    % Find change of sign.
    twosqrt = sqrt(2); 
    a = x; fa = fx; b = x; fb = fx;
    


    while (fa > 0) == (fb > 0)
        dx = twosqrt*dx;
        a = x - dx;  [fa, una, duga, dxpa, dxma, dypa, dyma] = Vdiff_solv(a,uT,uTinv,uxp,uyp,uxm,uym,g,params,masks);

        if ~isfinite(fa) || ~isreal(fa) || ~isfinite(a)
            b = NaN;
            return
        end

        if (fa > 0) ~= (fb > 0) % check for different sign
            % Before we exit the while loop, print out the latest interval
 
            break
        end
        
        b = x + dx;  [fb, un, dug, dxp, dxm, dyp, dym] = Vdiff_solv(b,uT,uTinv,uxp,uyp,uxm,uym,g,params,masks);
 
   

    end % while
end % if (numel(x) == 2)

fc = fb; unc = un; dugc = dug; dxpc = dxp; dxmc = dxm; dypc = dyp; dymc = dym;
% Main loop, exit from middle of the loop
while fb ~= 0 && a ~= b
    % Insure that b is the best result so far, a is the previous
    % value of b, and c is on the opposite side of the zero from b.
    if (fb > 0) == (fc > 0)
        c = a;  fc = fa;
        d = b - a;  e = d;
        unc = una; dugc = duga; dxpc = dxpa; dxmc = dxma; dypc = dypa; dymc = dyma;
    end
    if abs(fc) < abs(fb)
        a = b;    b = c;    c = a;
        fa = fb;  fb = fc;  fc = fa;
        una = un; duga = dug; dxpa = dxp; dxma = dxm; dypa = dyp; dyma = dym;
        un = unc; dug = dugc; dxp = dxpc; dxm = dxmc; dyp = dypc; dym = dymc;
        unc = una; dugc = duga; dxpc = dxpa; dxmc = dxma; dypc = dypa; dymc = dyma;
    end
    
    % Convergence test and possible exit
    m = 0.5*(c - b);
    toler = 2.0*tol*max(abs(b),1.0);
    if (abs(m) <= toler) || (fb == 0.0) 
        break
    end
   
    
    % Choose bisection or interpolation
    if (abs(e) < toler) || (abs(fa) <= abs(fb))
        % Bisection
        d = m;  e = m;
    else
        % Interpolation
        s = fb/fa;
        if (a == c)
            % Linear interpolation
            p = 2.0*m*s;
            q = 1.0 - s;
        else
            % Inverse quadratic interpolation
            q = fa/fc;
            r = fb/fc;
            p = s*(2.0*m*q*(q - r) - (b - a)*(r - 1.0));
            q = (q - 1.0)*(r - 1.0)*(s - 1.0);
        end
        if p > 0
            q = -q;
        else
            p = -p;
        end
        % Is interpolated point acceptable
        if (2.0*p < 3.0*m*q - abs(toler*q)) && (p < abs(0.5*e*q))
            e = d;  d = p/q;
        else
            d = m;  e = m;
        end
    end % Interpolation
    
    % Next point
    a = b;
    fa = fb;
    una = un; duga = dug; dxpa = dxp; dxma = dxm; dypa = dyp; dyma = dym;
    if abs(d) > toler
        b = b + d;
    elseif b > c
        b = b - toler;
    else
        b = b + toler;
    end
    [fb, un, dug, dxp, dxm, dyp, dym] = Vdiff_solv(b,uT,uTinv,uxp,uyp,uxm,uym,g,params,masks);
end % Main loop

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

function [V,D] = eignvec3x3S(A)
a1 = A(1);
a2 = A(2);
a3 = A(3);
a4 = A(4);
a5 = A(5);
a6 = A(6);

A1 = [a1 a2 a3;
      a2 a4 a5;
      a3 a5 a6];

p1 = a2^2 + a3^2 + a5^2;
q = (a1 + a4 + a6)/3;

a1 = a1 - q;
a4 = a4 - q;
a6 = a6 - q;

p2 = a1^2 + a4^2 + a6^2 + 2*p1;
p = sqrt(p2 / 6);

% Determinant
d = a1*(a5*a5 - a4*a6) + a2*(a2*a6-a5*a3) + a3*(a4*a3 - a2*a5);

r = -d/(2*p*p*p);

% In exact arithmetic for a symmetric matrix  -1 <= r <= 1
% but computation error can leave it slightly outside this range.
if (r <= -1)
    phi = pi / 3;
elseif (r >= 1)
    phi = 0;
else
    phi = acos(r) / 3;
end

% the eigenvalues satisfy eig3 <= eig2 <= eig1
eig1 = q + 2*p*cos(phi);
eig3 = q + 2*p*cos(phi + (2*pi/3));
eig2 = 3*q - eig1 - eig3;

if abs(eig1 - eig3) < 1E-14 || abs(eig1 - eig2) < 1E-14 || abs(eig2 - eig3) < 1E-14
    A = [A(1) A(2) A(3);
         A(2) A(4) A(5);
         A(3) A(5) A(6)];
    [V,D] = eig(A);
    D(1,1) = 1/realsqrt(D(1,1));
    D(2,2) = 1/realsqrt(D(2,2));
    D(3,3) = 1/realsqrt(D(3,3));
else
    A2 = A1;
    A3 = A1(:,1);
    V = A1;
    
    A1(1,1) = A1(1,1) - eig1;
    A1(2,2) = A1(2,2) - eig1;
    A1(3,3) = A1(3,3) - eig1;
    
    A2(1,1) = A2(1,1) - eig2;
    A2(2,2) = A2(2,2) - eig2;
    A2(3,3) = A2(3,3) - eig2;
    
    A3(1) = A3(1) - eig3;
    
    V(:,1) = A2*A3;
    % V(:,2) = A1*A3;
    V(:,3) = A1*A2(:,1);
    V(1,2) = V(2,1)*V(3,3) - V(3,1)*V(2,3);
    V(2,2) = V(3,1)*V(1,3) - V(1,1)*V(3,3);
    V(3,2) = V(1,1)*V(2,3) - V(2,1)*V(1,3);
    
    V(:,1) = V(:,1)/realsqrt(V(1,1)^2 + V(2,1)^2 + V(3,1)^2);
    V(:,2) = V(:,2)/realsqrt(V(1,2)^2 + V(2,2)^2 + V(3,2)^2);
    V(:,3) = V(:,3)/realsqrt(V(1,3)^2 + V(2,3)^2 + V(3,3)^2);
    
    D = [1/realsqrt(eig1), 1/realsqrt(eig2), 1/realsqrt(eig3)];
end
end