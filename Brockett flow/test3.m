n = 20;
rng(10)
Q = rand(n);
Q = (Q + Q')/2;
Q = 4*(rand(n,1)-0.5);
P = orth(rand(n));
Q = P'*diag(Q)*P;
[V,D] = eig(Q);
[DD,I] = sort(diag(D));
V = V(:,I);


N = length(Q):-1:1;
N = diag(N);
theta = eye(n)
f = @(th) -th*(N*th'*Q*th - th'*Q*th*N);
F = @(th) trace(Q*th*N*th');

dt = 0.01;

M = 120
Fdiff = zeros(M,1);

for i = 1:M
    w = zeros(n);
    theta0 = theta;
    for k = 1:n
        for l = k+1:n
            ekl = zeros(n);
            ekl(k,l) = 1;
            ekl(l,k) = -1;
            G = @(a) a + dt*(F(expm(w + a*ekl)*theta0) - F(theta))/a;
            alpha = fzero(G,0.01);
            w = w + alpha*ekl;
            theta = expm(w)*theta0;
        end
    end
    if ~mod(i,10)
        QQ = theta'*Q*theta;
%         F(theta) - F(V)
% diag(QQ)
% diag(D)
%         norm(diag(QQ) - diag(D))
    end
        Fdiff(i) = F(theta) - F(V);
%     Fdiff(i) = norm(diag(theta'*Q*theta) - DD);
    norm(diag(theta'*Q*theta) - DD)
    Fdiff(i)
%     dt = dt*i/(i+1);
end

figure
semilogy(dt*(1:M),Fdiff,'linewidth',1.5)
hold on
semilogy(dt*(1:M),1./(dt*(1:M)).^2,'linewidth',1.5)
semilogy(dt*(1:M),1./(dt*(1:M)),'linewidth',1.5)
semilogy(dt*(1:M),exp(-0.25*dt*(1:M)),'linewidth',1.5)
legend('DG','1/n^2','1/n','e^{-0.025t}','location','northeast')
xlabel('t')
ylabel('Optimality error')
set(gcf,'paperposition',[0 0 5 4.2])
set(gca,'linewidth',1);
filename = ['BrockettFlowConvExp'];
print('-depsc','-r300',[filename '.eps'])
