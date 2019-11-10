function g = itab(u,V,alpha)
n = length(u);
k = 0;
wj = u;
V_old = V(u);
g = zeros(n);
for i = 1:n
    for j = i:n
        k = k + 1;
        ekl = zeros(n);
        ekl(i,j) = 1;
        ekl(j,i) = 1;
        wj = wj + alpha(k)*ekl;
        V_new = V(wj);       
        g = g + (V_new - V_old)/alpha(k)*ekl;
        V_old = V_new;
    end
end
end

