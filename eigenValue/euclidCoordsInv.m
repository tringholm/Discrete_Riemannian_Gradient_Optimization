function phi = euclidCoordsInv(x)
n = length(x);
phi = zeros(n-1,1);

for i = 1:n-2
    phi(i) = acos(x(i)/sqrt(x(i:end)'*x(i:end)));
end
if x(end) < 0
    phi(end) = 2*pi - acos(x(n-1)/sqrt(x(n)^2 + x(n-1)^2));
else
    phi(end) = acos(x(n-1)/sqrt(x(n)^2 + x(n-1)^2));
end
end