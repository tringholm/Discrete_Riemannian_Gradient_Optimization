function v = getTri(A)
N = length(A);
v = zeros(N*(N+1)/2,1);
k = 1;
for i = 1:N
    v(k:k+N-i) = A(i:N,i);
    k = k + N + 1 - i;
end
end