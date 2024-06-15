function [nab4,P,D]=disc4(cfs)
a=cfs(1);
b=cfs(2);
c=cfs(3);
d=cfs(4);
e=cfs(5);
nab4=256.*a.^3.*e.^3-192.*a.^2.*b.*d.*e.^2-128.*a.^2.*c.^2.*e.^2+144.*a.^2.*c.*d.^2.*e-27.*a.^2.*d.^4+144.*a.*b.^2.*c.*e.^2-6.*a.*b.^2.*d.^2.*e ...
-80.*a.*b.*c.^2.*d.*e+18.*a.*b.*c.*d.^3+16.*a.*c.^4.*e-4.*a.*c.^3.*d.^2-27.*b.^4.*e.^2+18.*b.^3.*c.*d.*e-4.*b.^3.*d.^3-4.*b.^2.*c.^3.*e ...
+b.^2.*c.^2.*d.^2;
nab4=simplify(nab4);
P=8*a.*c-3*b.^2;
D=64*a.^3.*e-16*a.^2.*c.^2+16*a.*b.^2.*c-16*a.^2.*b.*d-3*b.^4;
end