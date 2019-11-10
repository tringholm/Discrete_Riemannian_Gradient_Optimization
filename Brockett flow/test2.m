n = 7;
H = rand(n);
H = (H + H')/2;
eig(H)

N = 1:length(H);
N = diag(N);

brack = @(x,y) x*y - y*x;
f = @(x) brack(x,brack(x,N));
F = @(x) trace((x*N - N*x)*(x*N - N*x)')

dt = 1E-2;

H0 = H;
M = 10000
for i = 1:M
    for k = 1:n
        for l = 1:k
            ekl = zeros(n);
            ekl(k,l) = 1;
            ekl(l,k) = 1;
            G = @(a) a + dt*(F(H + a*ekl) - F(H))/a;
            alpha = fzero(G,0.01);
            if k == l
               alpha 
            end
            H = H + alpha*ekl;
        end
    end
    if ~mod(i,100)
        H
        F(H)
    end
end