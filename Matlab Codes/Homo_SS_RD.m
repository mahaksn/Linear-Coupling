%% Homogeneous Steady state

function [u_star,v_star]=Homo_SS_RD(f,g,a,b)

syms u v

ff=f(u,v,a,b);
gg=g(u,v,a,b);

vv=solve(gg==0,v);
ff1=subs(ff,v,vv);
ff2=solve(ff1==0,u);

u_star=ff2(1);
v_star=subs(vv,u,u_star);

end
