% Denoise InSAR image 

function [u,Vhist] = InSAR_DG_NoSmooth_timing(g,p,q,gamma,tol,xtol,dt,Vmin,Vminfact)
u = g;
% u = mod(g + randn(size(g))+pi,2*pi)-pi;

Tmax = 20000;
residual = inf;

Vhist = zeros(Tmax,1);
Vhist(1) = Veval(u,g,p,q,gamma);
VV = Vhist(1);
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

% doit = 1;
for k = 1:Tmax
    if ~mod(k,200)
        dt = dt/2
    end
   if  Vhist(k) < Vmin*Vminfact || residual < 0
       k
      break 
   end
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
%            solveFxn = @(alp) alp + dt*Vdiff(alp,uo,uxp,uyp,uxm,uym,g(i,j),p,q,gamma/q,i,j,M,N,Vo)/alp;
%            [alpha,fval] = fzero_noExtra(solveFxn,0.01,xtol);%,optimset('Display','iter'));
           [alpha,fval] = fzero_noExtra_spec([(-2*pi -uo), (2*pi - uo )],uo,uxp,uyp,uxm,uym,g(i,j),p,gamma/q,dt,i,j,M,N,Vo,xtol);
           VV = VV - alpha^2/dt;
%            [alpha,fval] = fzero_noExtra_spec(0.01,uo,uxp,uyp,uxm,uym,g(i,j),p,q,gamma/q,dt,i,j,M,N,Vo,xtol);
%            if alpha < -pi -uo 
%                alpha + pi + uo
%            end
%            if alpha > pi - uo
%                alpha - pi + uo
%            end
%             lval = Vdiff_solv((-3*pi -uo),uo,uxp,uyp,uxm,uym,g(i,j),p,q,gamma/q,dt,i,j,M,N,Vo);
%             rval = Vdiff_solv((3*pi - uo ),uo,uxp,uyp,uxm,uym,g(i,j),p,q,gamma/q,dt,i,j,M,N,Vo);
% 
%            al = -3*pi - uo:0.001:pi - uo;
%            y = zeros(size(al));
% 
%            for k = 1:length(y)
%                y(k) = Vdiff_solv(al(k),uo,uxp,uyp,uxm,uym,g(i,j),p,q,gamma/q,dt,i,j,M,N,Vo);
%            end
%            plot(al,y)
%            if lval*rval > 0
%                lval
%                rval
%            end
%            if fval> 1E-2 && alpha > 100*xtol
%                fval
%                alpha
%            end
           u_new = mod(u(i,j) + alpha + pi,2*pi) - pi;

%            u_new = rem(uo+alpha,pi);

           u(i,j) = u_new;
           dp(i,j) = d_ang(u_new,g(i,j))^p;
           if j > 1
               dqx(i,j) = (pi - abs(abs(u_new-uxm) -pi));
           end
           if j < N
               dqx(i,j+1) = (pi - abs(abs(u_new-uxp) -pi));
           end
           if i > 1
               dqy(i,j) = (pi - abs(abs(u_new-uym) -pi));
           end
           if i < M
               dqy(i+1,j) = (pi - abs(abs(u_new-uyp) -pi));
           end
       end
   end
   Vhist(k+1) = Veval(u,g,p,q,gamma);
%    Vhist(k+1) = VV;
%    VV - Vhist(k+1)
%    figure(imfig); imagesc(u); title('Intermediate plot');
%    pause(0.01)
    residual = (Vhist(k) - Vhist(k+1))/Vhist(1);
    if ~mod(k, 10)
        n = 1:10:k+1;
        figure(imfig); loglog(n,(Vhist(1:10:k+1) - Vmin)/(Vhist(1) - Vmin));
        hold on
        plot(n,1./n.^0.5)
        plot(n,1./n)
        plot(n,1./n.^2)
        drawnow()
        residual
        hold off
    end
end
residual
Vhist = Vhist(~~Vhist);
end


function out = Vdiff_solv(alp,u,uxp,uyp,uxm,uym,g,p,gamma,dt,i,j,M,N,Vo)
un = u + alp + pi;
if un > 2*pi
    un = mod(un,2*pi) - pi;
elseif un < 0
    un = mod(un,2*pi) - pi;
else
    un = un - pi;
end
if p == 1
    out = (pi - abs(abs(un-g) - pi));
elseif p == 2
    out = 0.5*(pi - abs(abs(un-g) - pi))^2;
end
out2 = (j<N)*(pi - abs(abs(un-uxp) -pi)) + (i<M)*(pi - abs(abs(un-uyp) -pi)) + (j>1)*(pi - abs(abs(un-uxm) -pi)) + (i>1)*(pi - abs(abs(un-uym) -pi));
out = out + gamma*out2 - Vo;
out = alp + dt*out/alp;
end



function out = d_ang(phi,theta)
% out = pi - abs(abs(phi-theta) -pi); 
out = pi - abs(abs(phi-theta) - pi);
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



function [b,fval] = fzero_noExtra_spec(x,u,uxp,uyp,uxm,uym,g,p,gamma,dt,i,j,M,N,Vo,tol)

% Initialization
fcount = 0;
iter = 0;
intervaliter = 0;


% Interval input
if (numel(x) == 2) 
 
    a = x(1);
    b = x(2);    
    fa = Vdiff_solv(a,u,uxp,uyp,uxm,uym,g,p,gamma,dt,i,j,M,N,Vo);
  
    fb = Vdiff_solv(b,u,uxp,uyp,uxm,uym,g,p,gamma,dt,i,j,M,N,Vo);
 
    fcount = fcount + 2;
    
    if ( fa == 0 )
        b = a;
        fval = fa;
        return
    elseif ( fb == 0)
        % b = b;

        fval = fb;
        return
    end
    
    % Starting guess scalar input
elseif (numel(x) == 1)

    
    fx = Vdiff_solv(x,u,uxp,uyp,uxm,uym,g,p,gamma,dt,i,j,M,N,Vo);

    fcount = fcount + 1;  
    if fx == 0
        b = x;
 
        fval = fx;
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
        intervaliter = intervaliter + 1;
        dx = twosqrt*dx;
        a = x - dx;  fa = Vdiff_solv(a,u,uxp,uyp,uxm,uym,g,p,gamma,dt,i,j,M,N,Vo);
        fcount = fcount + 1;
        if ~isfinite(fa) || ~isreal(fa) || ~isfinite(a)
            b = NaN; fval = NaN;
            return
        end

        if (fa > 0) ~= (fb > 0) % check for different sign
            % Before we exit the while loop, print out the latest interval
 
            break
        end
        
        b = x + dx;  fb = Vdiff_solv(b,u,uxp,uyp,uxm,uym,g,p,gamma,dt,i,j,M,N,Vo);
 
        fcount = fcount + 1;        
   

    end % while
end % if (numel(x) == 2)

fc = fb;
% Main loop, exit from middle of the loop
while fb ~= 0 && a ~= b
    % Insure that b is the best result so far, a is the previous
    % value of b, and c is on the opposite side of the zero from b.
    if (fb > 0) == (fc > 0)
        c = a;  fc = fa;
        d = b - a;  e = d;
    end
    if abs(fc) < abs(fb)
        a = b;    b = c;    c = a;
        fa = fb;  fb = fc;  fc = fa;
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
            pp = 2.0*m*s;
            qq = 1.0 - s;
        else
            % Inverse quadratic interpolation
            qq = fa/fc;
            r = fb/fc;
            pp = s*(2.0*m*qq*(qq - r) - (b - a)*(r - 1.0));
            qq = (qq - 1.0)*(r - 1.0)*(s - 1.0);
        end
        if pp > 0
            qq = -qq;
        else
            pp = -pp;
        end
        % Is interpolated point acceptable
        if (2.0*pp < 3.0*m*qq - abs(toler*qq)) && (pp < abs(0.5*e*qq))
            e = d;  d = pp/qq;
        else
            d = m;  e = m;
        end
    end % Interpolation
    
    % Next point
    a = b;
    fa = fb;
    if abs(d) > toler
        b = b + d;
    elseif b > c
        b = b - toler;
    else
        b = b + toler;
    end
    fb = Vdiff_solv(b,u,uxp,uyp,uxm,uym,g,p,gamma,dt,i,j,M,N,Vo);
    fcount = fcount + 1;
    iter = iter + 1;
end % Main loop

fval = fb; % b is the best value


end

