function del3=disc3(cfs)
a=cfs(1);
b=cfs(2);
c=cfs(3);
d=cfs(4);
del3=b^2*c^2-4*a*c^3-4*b^3*d-27*a^2*d^2+18*a*b*c*d;
del3=simplify(del3);