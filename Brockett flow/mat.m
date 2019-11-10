function w = mat(alpha,n)
w = zeros(n);
m = 0;
for k = 1:n
    for l = k:n
        m = m + 1;
        w(k,l) = alpha(m);
        w(l,k) = alpha(m);
    end
end