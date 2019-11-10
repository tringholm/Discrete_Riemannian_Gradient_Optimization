function out = Vdiff_solv_func(alp,u,uxp,uyp,uxm,uym,g,p,q,gamma,dt,i,j,M,N,Vo)
% un = rem(u+alp,pi);
un = mod(u + alp + pi,2*pi) - pi;
out = 1/p*(pi - abs(abs(un-g) -pi))^p;
% out2 = (j<N)*d_ang(un,uxp) + (j>1)*d_ang(uxm,un) + (i<M)*d_ang(un,uyp) + (i>1)*d_ang(uym,un);
out2 = (j<N)*(pi - abs(abs(un-uxp) -pi)) + (j>1)*(pi - abs(abs(un-uxm) -pi)) + (i<M)*(pi - abs(abs(un-uyp) -pi)) + (i>1)*(pi - abs(abs(un-uym) -pi));
out = out + gamma*out2 - Vo;
out = alp + dt*out/alp;
end

