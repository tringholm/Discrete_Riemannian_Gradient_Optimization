close all
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
theta = eye(n);
thetaexp = theta;
f = @(th) -th*(N*th'*Q*th - th'*Q*th*N);
F = @(th) trace(Q*th*N*th');
I = eye(length(N));
cay = @(A) (I+A)\(I-A);

dt = 0.1;

M = 1000
Fdiff = zeros(M,1);

eigens = zeros(n,M);
for i = 1:M
    w = zeros(n);
%     wexp = zeros(n);
    theta0 = theta;
%     theta0exp = thetaexp;
    for k = 1:n
        for l = k+1:n
            ekl = zeros(n);
            ekl(k,l) = 1;
            ekl(l,k) = -1;
            
%             ek = zeros(n,1);
%             el = zeros(n,1);
%             ek(k) = 1;
%             el(l) = 1;
%             ekl = [ek el]*[el';-ek']
            
            
            
            G = @(a) a + dt*(F(cay(w + a*ekl)*theta0) - F(theta))/a;
%             H = @(a) a + dt*(F(expm(wexp + a*ekl)*theta0exp) - F(thetaexp))/a;
            alpha = fzero(G,0.01);
%             alphaexp = fzero(H,0.01);
% phi1 = cay(w);
% phi2 = cay(w + alpha*ekl);
% phi2 - phi1
            w = w + alpha*ekl;
%             wexp = wexp + alphaexp*ekl;
            theta = cay(w)*theta0;
%             thetaexp = expm(wexp)*theta0exp;
        end
    end
    eigens(:,i) = diag(theta'*Q*theta);
    if ~mod(i,10)
        QQ = theta'*Q*theta;
        F(theta) - F(V)
%         F(thetaexp) - F(V)
% diag(QQ)
% diag(D)
%         norm(diag(QQ) - diag(D))
    end
%     dt = dt*i/(i+1);
    Fdiff(i) = F(theta) - F(V);
    Fdiff(i) = norm(diag(theta'*Q*theta) - DD);
    norm(diag(theta'*Q*theta) - DD)
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
filename = ['BrockettFlowConvEigen'];
print('-depsc','-r300',[filename '.eps'])

figure
for i = 1:n
   plot(dt*(1:M),eigens(i,:),'linewidth',1.5)
   hold on
end

ax = gca;
ax.ColorOrderIndex = 1;


for i = 1:n
   plot(dt*(1:10:(M+1)),DD(i)*ones(length((1:10:(M+1))),1),':','linewidth',2.5)
   hold on
end

axis([0 dt*M -2 2])

xlabel('t')
ylabel('Eigenvalues ')
set(gcf,'paperposition',[0 0 5 4.2])
set(gca,'linewidth',1);
filename = ['BrockettFlow'];
print('-depsc','-r300',[filename '.eps'])