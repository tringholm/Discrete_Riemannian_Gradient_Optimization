function out = Vdiff(alp,u,uxp,uyp,uxm,uym,g,p,q,gamma,i,j,M,N,Vo)
un = mod(u + alp + pi,2*pi) - pi;
out = 1/p*d_ang_func(un,g)^p;
if j < N
    out = out + (j < N)*gamma*d_ang_func(un,uxp)^q ;
end
if j > 1
    out = out + gamma*d_ang_func(uxm,un)^q;
end
if i < M
    out = out + gamma*d_ang_func(un,uyp)^q;
end
if i > 1
    out = out + gamma*d_ang_func(uym,un)^q;
end
out = out - Vo;
end

