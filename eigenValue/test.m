clear all
close all
rng(20);

N = 400;

A = rand(N);
A = (A+A')/2;

dt = 1;
tol = 1E-10;
u = 1 - 2* rand(N,1);
u = u/norm(u);

phi = euclidCoordsInv(u);

% tic;
% v1=eigenValueSphere2(A,phi,tol,dt);
% toc
tic;
v1=eigenValueSphere3(A,phi,tol,dt);
toc
v1(end)
tic;
v2=eigenValueSphere3_mex(A,phi,tol,dt);
toc
v2(end)
% v2 - v1