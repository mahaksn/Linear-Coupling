function [Jac,stable]=stability(f,g,u,v,u0,v0)
Jac = jacobian([f,g],[u,v]);
Jac = double(subs(Jac,[u,v],[u0,v0]));

CPJac=charpoly(Jac);
[rhCPJac,stable]=rhSCD_sym(CPJac);

if stable==1
    fprintf('\nThe system is stable without diffusion.\n\n');
elseif stable==0
    fprintf('\nThe system is unstable without diffusion.\n\n');
end
end
