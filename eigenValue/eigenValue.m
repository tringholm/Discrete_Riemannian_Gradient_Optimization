% Compute the lowest eigenvalue/vector of a symmetric matrix by minimizing the
% Rayleigh quotient over S^{n-1}. 
clear all
close all
rng(20);

N = 20;

A = rand(N);
A = (A+A')/2;

Veval = @(z) z'*A*z;
getBasis = @(z) null(z');
basePoint = @(z) z;
retraction = @(p,v_p) (p+v_p)/norm(p+v_p);
retractionInv = @(p,q) q/(p'*q) - p;
% Vdiff = @(c_p,p,v_p,alp,ret,E) Veval(ret(c_p,v_p + alp*E)) - Veval(p);
Vdiff = @(c_p,p,v_p,alp,ret,E) (ret(c_p,v_p + alp*E)-p)'*A*(ret(c_p,v_p + alp*E)+p);


dt = 0.25;
tol = 1E-12;
u = rand(N,1);
u = u/norm(u);

Tmax = 1000;
residual = inf;

Vhist = zeros(Tmax,1);
Vhist(1) = Veval(u);
tic
for i = 1:Tmax
   if  residual < tol
      break 
   end
   
   c = basePoint(u);
   E = getBasis(c);
   v = u;
   w = retractionInv(c,v);
   
   N = size(E,2);
   if mod(i,2)
       elems = 1:N;
   else
       elems = N:-1:1;
   end
   
   for j = elems
       Ej = E(:,j);
       solveFxn = @(alp) alp + dt*Vdiff(c,v,w,alp,retraction,Ej)/alp;
       alpha = fzero(solveFxn,0.01);
       w = w + alpha*Ej;
       v = retraction(c,w);
   end
   u = v;
   Vhist(i+1) = Veval(u);
%    Vhist(i+1)
   residual = (Vhist(i) - Vhist(i+1))/Vhist(1);
end
toc
residual
Vhist = Vhist(~~Vhist);
i
[V,D]= eig(A);
min(norm(V(:,1)-u),norm(V(:,1)+u))
D(1,1) - Vhist(end)