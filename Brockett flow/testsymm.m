n = 5;
rng(15)
Q = rand(n);
Q = (Q + Q')/2;
[V,D] = eig(Q);


N = length(Q):-1:1;
N = diag(N);
H = Q;

F = @(th) trace(th*N);
brack = @(x,y) x*y -y*x;
R = @(u,x) brack(u,brack(u,x));


dt = 0.1;

M = 100000;
nn = n*(n+1)/2;
for i = 1:M
    G = @(a) getTri(mat(a,n)) + dt*getTri(R(H,itab(H,F,a)));
    options = optimoptions('fsolve','Display','off');
    alpha = fsolve(G,zeros(nn,1) + 0.01);
    ss = dt*R(H,N);
    H - ss
    H = H + mat(alpha,n)
    alpha;
    if ~mod(i,10)
        QQ = H'*Q*H;
        F(H) - F(D)
%         F(thetaexp) - F(V)
% diag(QQ)
% diag(D)
%         norm(diag(QQ) - diag(D))
    end
%     dt = dt*i/(i+1);
end