function [V,D] = eignvec3x3(A)
a1 = A(1,1);
a2 = A(2,1);
a3 = A(3,1);
a4 = A(2,2);
a5 = A(3,2);
a6 = A(3,3);

p1 = a2^2 + a3^2 + a5^2;
q = (a1 + a4 + a6)/3;

a1 = a1 - q;
a4 = a4 - q;
a6 = a6 - q;

p2 = a1^2 + a4^2 + a6^2 + 2*p1;
p = sqrt(p2 / 6);

% Determinant
d = a1*(a5*a5 - a4*a6) + a2*(a2*a6-a5*a3) + a3*(a4*a3 - a2*a5);

r = -d/(2*p*p*p);

% In exact arithmetic for a symmetric matrix  -1 <= r <= 1
% but computation error can leave it slightly outside this range.
if (r < -1)
    phi = pi / 3;
elseif (r > 1)
    phi = 0;
else
    phi = acos(r) / 3;
end

% the eigenvalues satisfy eig3 <= eig2 <= eig1
eig1 = q + 2*p*cos(phi);
eig3 = q + 2*p*cos(phi + (2*pi/3));
eig2 = 3*q - eig1 - eig3;
D = [1/realsqrt(eig1), 1/realsqrt(eig2), 1/realsqrt(eig3)];

A1 = A;
A2 = A;
A3 = A(:,1);
V = A;

A1(1,1) = A1(1,1) - eig1;
A1(2,2) = A1(2,2) - eig1;
A1(3,3) = A1(3,3) - eig1;

A2(1,1) = A2(1,1) - eig2;
A2(2,2) = A2(2,2) - eig2;
A2(3,3) = A2(3,3) - eig2;

A3(1) = A3(1) - eig3;

V(:,1) = A2*A3;
V(:,2) = A1*A3;
V(:,3) = A1*A2(:,1);

V(:,1) = V(:,1)/realsqrt(V(1,1)^2 + V(2,1)^2 + V(3,1)^2);
V(:,2) = V(:,2)/realsqrt(V(1,2)^2 + V(2,2)^2 + V(3,2)^2);
V(:,3) = V(:,3)/realsqrt(V(1,3)^2 + V(2,3)^2 + V(3,3)^2);
end