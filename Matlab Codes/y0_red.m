function sss=y0_red(y0h,h,k,eta,k0,eta0)
y0k=subs(y0h,h,k^2);
e=subs(y0k,eta,eta0);
sss=subs(e,k,k0);
end
