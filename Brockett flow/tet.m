H = rand(7)
H = (H + H')/2;
J = eig(H)

N = 1:length(H);
N = diag(N);

brack = @(x,y) x*y - y*x;
f = @(x) brack(x,brack(x,N));

dt = 1E-5;

H0 = H;
M = 10000000
for i = 1:M
    H = H - dt*f(H);
    if ~mod(i,100)
        H
        trace((H*N - N*H)*(H*N - N*H)')
    end
end