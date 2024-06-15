function [RHv,k00,eta00]=ycoeff(y,k0,eta0)
[k00,eta00]=meshgrid(k0,eta0);
RH=vectorize(y);
RH=str2func(['@(k,eta)' RH])
RHv=RH(k00,eta00);
end
