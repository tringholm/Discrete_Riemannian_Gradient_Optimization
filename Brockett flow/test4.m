n = 7;
rng(12)
Q = rand(n);
Q = (Q + Q')/2;
J = eig(Q)

N = 1:length(Q);
N = diag(N);
theta = eye(n)
f = @(th) -th*(N*th'*Q*th - th'*Q*th*N);
F = @(th) trace(Q*th*N*th');

dt = 1E-4;

M = 10000000
for i = 1:M
    theta = theta - dt*f(theta);
    if ~mod(i,100)
        theta'*Q*theta
        F(theta)
    end
end