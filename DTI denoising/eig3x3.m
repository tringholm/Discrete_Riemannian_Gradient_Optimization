function [eig1, eig2, eig3] = eig3x3(a1,a2,a3,a4,a5,a6)
p1 = a2^2 + a3^2 + a5^2;
q = (a1 + a4 + a6)/3;

a1 = a1 - q;
a4 = a4 - q;
a6 = a6 - q;

p2 = a1^2 + a4^2 + a6^2 + 2*p1;
p = 2*sqrt(p2 / 6);

% Determinant
d = a1*(a5*a5-a4*a6) + a2*(a2*a6-a5*a3) + a3*(a4*a3-a2*a5);

r = -4*d/(p*p*p);

% In exact arithmetic for a symmetric matrix  -1 <= r <= 1
% but computation error can leave it slightly outside this range.
if (r < -1)
    phi = pi/3;
elseif (r > 1)
    phi = 0;
else
    phi = acos(r)/3;
end

% the eigenvalues satisfy eig3 <= eig2 <= eig1
eig1 = q + p*cos(phi);
eig3 = q + p*cos(phi + (2*pi/3));
eig2 = 3*q - eig1 - eig3;
% out = [eig1, eig2, eig3];
end