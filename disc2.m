function del2=disc2(cfs)
a=cfs(1);
b=cfs(2);
c=cfs(3);
del2=b^2-4*a*c;
del2=simplify(del2);