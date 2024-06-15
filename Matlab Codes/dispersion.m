function Mr=dispersion(CP,k,eta,k0,eta0)
[k00,eta00]=meshgrid(k0,eta0);
for j=1:length(eta00)
    for i=1:length(k00)
    CP1=subs(CP,[k,eta],[k00(j,i),eta00(j,i)]);
    EG1=solve(CP1);
    rp=real(EG1);
    % img=imag(EG1);
    Mr(j,i)=max(rp);
    % Mi(j,i)=max(img);
    end
end
Mr=double(Mr);
% Mi=double(Mi);
end
