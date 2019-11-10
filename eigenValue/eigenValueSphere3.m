% Compute the lowest eigenvalue/vector of a symmetric matrix by minimizing the
% Rayleigh quotient over S^{n-1}.

function Vhist = eigenValueSphere3(A,phi,tol,dt)
%#codegen
N = length(phi)+1;
Tmax = 2000;
residual = inf;
tolNewton =  1E-13;

dA = diag(A);
AU = triu(A,1);
Vhist = zeros(Tmax,1);
u = euclidCoords(phi);
rowsums = u.*(A*u);
Vhist(1) = sum(rowsums);
% tic
sums = zeros(1,5);
sums(1) = A(1,1)*u(1)^2;                          % i = k, j = k
sums(2) = 0;                                      % i < k, j = k
sums(3) = 2*(rowsums(1) - sums(1));               % i = k, j > k
sums(4) = 0;                                      % i < k, j > k
sums(5) = Vhist(1) - sums(3) - sums(1);           % i > k, j > k
trisum = 2*u.*(AU*u);

newtoniters = 0;
t = 0;
for i = 1:Tmax
    if  residual < tol
        break
    end
    sum3fact = 1;
    ufact = 1;
    for j = 1:N-1
    
        
        alpha = 0.1;
        dta = dt/alpha;
        cp = 1/cos(phi(j)); sp = 1/sin(phi(j));
        cpa = cos(phi(j)+alpha); spa = sin(phi(j)+alpha);
        cc = cpa*cp; ss = spa*sp; sc = spa*cp; cs = cpa*sp;
        
        Vd = sums(1)*(cc^2-1) + sums(2)*(cc-1) + sums(3)*(cc*ss - 1) + sums(4)*(ss-1) + sums(5)*(ss^2-1);
        V = alpha + dta*Vd;
        V_old = inf;
        cnt = 0;
        while abs(V) > tolNewton && cnt < 10 && V_old - V > 0
            dV = sums(1)*(-2*sc*cc) + sums(2)*(-sc) + sums(3)*(cc*cs - sc*ss) + sums(4)*cs + sums(5)*2*cs*ss;
            dV = 1+ dta*dV - dta*Vd/alpha;
            alpha = alpha - V/dV;
            dta = dt/alpha;
            cpa = cos(phi(j)+alpha); spa = sin(phi(j)+alpha);
            cc = cpa*cp; ss = spa*sp; sc = spa*cp; cs = cpa*sp;
            Vd = sums(1)*(cc^2-1) + sums(2)*(cc-1) + sums(3)*(cc*ss-1) + sums(4)*(ss-1) + sums(5)*(ss^2-1);
            V_old = V;
            V = alpha + dta*Vd;
            cnt = cnt + 1;
        end
        
        
        newtoniters = newtoniters + cnt;

        u(j) = u(j)*cc*ufact;
        ufact = ufact*ss;
        un = u(j+1)*ufact;
        sum3fact = sum3fact*ss^2;
        
        sums(1) = un*A(j+1,j+1)*un;              % i = k, j = k
        sums(2) = 0;

        for k = 1:j
            sums(2) = sums(2) + A(k,j+1)*u(k);
        end
        sums(2) = 2*un*sums(2);
        sums(4) = cc*ss*sums(3) + ss*sums(4) - sums(2);           % i < k, j > k
        sums(3) = trisum(j+1)*sum3fact;
        sums(5) = ss*ss*sums(5) - sums(3) - sums(1);           % i > k, j > k
        phi(j) = phi(j) + alpha;
    end
    u(end) = u(end)*ufact;
    trisum = 2*u.*(AU*u);

    dsum = sum(u.*dA.*u);
    Vhist(i+1) = sum(trisum)+dsum;
    sums(1) = A(1,1)*u(1)^2;                          % i = k, j = k
    sums(2) = 0;                                      % i < k, j = k
    sums(3) = trisum(1);                              % i = k, j > k
    sums(4) = 0;                                      % i < k, j > k
    sums(5) = Vhist(i+1) - sums(3) - sums(1);         % i > k, j > k
    residual = (Vhist(i) - Vhist(i+1))/abs(Vhist(1));
end
% toc
Vhist = Vhist(~~Vhist);
phi
end

function x = euclidCoords(phi)
n = length(phi);
x = zeros(n+1,1);
x(1) = cos(phi(1));
sinefact = sin(phi(1));
for i = 2:n
    x(i) = cos(phi(i))*sinefact;
    sinefact = sinefact*sin(phi(i));
end
x(end) = sinefact;
end