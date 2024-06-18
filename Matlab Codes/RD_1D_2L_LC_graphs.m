clc; clear; close all;
tic
% U-layer_1; L-layer_2;
% u,v - concentrations;

dim='1D';
layer='2L';
suffix='LC_graphs';
pdfname='RD_1D_2L_LC_graphs';

% format short
% format bank 

%% pre-setting
showanimation=1;
makegif=1;
drawperframe=200;
pattern=1; % Pattern fully developed at..

Schnakenberg=1;
CDIMA=0;
Brusselator=0;

case1=0;
case2=1;

%% set the parameters which govern equation
if Schnakenberg
RD='Schnak';
f = @(u,v,a,b) a - u + (u.^2).*v;
g = @(u,v,a,b) b - (u.^2).*v;
if case1
    RD=[RD,'_','case1'];
    a1=0.05; b1=1.4; 
    aU=a1; bU=b1; aL=a1; bL=b1;
    DuU=0.00025; DvU=0.0125;
    DuL=0.01; DvL=0.2;
    eta_end=10;
    kplus=60;
elseif case2
    RD=[RD,'_','case2'];
    a=sym(0.05); b=sym(1.4);
    aU=a; bU=b; aL=a; bL=b;
    DuU=sym(1); DvU=sym(50);
    DuL=sym(40); DvL=sym(800);
    L=130;
    eta_end=10;
    kplus=1;
end
end

%% saving

[filename,folder] = savingcode(RD,dim,suffix,layer);
prefix = strcat(filename);
prefix = strcat(folder, prefix);

%% saving code
options = struct('format','pdf', ...
    'outputDir',prefix, ...
    'codeToEvaluate','0;');
publish(pdfname,options);
diary([prefix,'.txt']);
fprintf('saving to %s\n',folder);

%% Steady states

[uU0,vU0]=Homo_SS_RD(f,g,aU,bU);
[uL0,vL0]=Homo_SS_RD(f,g,aL,bL);
uU0=double(uU0); vU0=double(vU0);
uL0=double(uL0); vL0=double(vL0);
fprintf(['\nInitial Steady States for\n' ...
    'layer_1 are: uU0=%.5f, vU0=%.5f\n'...
    'layer_2 are: uL0=%.5f, vL0=%.5f\n'], ...
    uU0,vU0,uL0,vL0);

%% Linear stability analysis
syms uUs vUs uLs vLs
syms eta y k h
sympref('FloatingPointOutput',true);

fUs = f(uUs,vUs,aU,bU);
gUs = g(uUs,vUs,aU,bU);
fLs = f(uLs,vLs,aL,bL);
gLs = g(uLs,vLs,aL,bL);

disp('Stability of layer_1')
[JacU,stableU]=stability(fUs,gUs,uUs,vUs,uU0,vU0);
if stableU==1
    fprintf('~~~~~> it is a stable system! <~~~~~\n')
else
    fprintf('~~~~~> it is an unstable system! <~~~~~\n')
end

disp('Stability of layer_2')
[JacL,stableL]=stability(fLs,gLs,uLs,vLs,uL0,vL0);
if stableL==1
    fprintf('~~~~~> it is a stable system! <~~~~~\n')
else
    fprintf('~~~~~> it is an unstable system! <~~~~~\n')
end

%% Coulped system with diffusion
DM = [DuU,0,0,0;
    0,DvU,0,0;
    0,0,DuL,0;
    0,0,0,DvL];
    
eta_uu=eta;
eta_vv=eta;
CUM = [eta_uu,0;0,eta_vv];
CLM = [eta_uu,0;0,eta_vv];
CM = [-CUM CUM; CLM -CLM];

if stableU==1 && stableL==1
    Jac = [JacU,zeros(2);zeros(2),JacL];
    scriptL = Jac - k^2*DM + CM;
    CPy = charpoly(scriptL,y);
    disp('The J_{\eta} matrix is:')
    disp(scriptL)
    disp('The Dispersion relation is:')
    disp(CPy)
    [cfsCPk,ysCPk]=coeffs(CPy,y);
    rhCPk=rhSCD_sym(cfsCPk);
    % disp('Routh Table for coupled system with diffusion')
    % disp(rhCPk)
end

%% Values of the Coefficients of the dispersion relation

disp('Coefficients of the dispersion relation:')
y3=cfsCPk(2)
y2=cfsCPk(3)
y1=cfsCPk(4)
y0=cfsCPk(5)
ybnum=simplify(y3*y2-y1);
yb=simplify(ybnum/y3)   % yb - ybeta
[ygN,ygD]=numden(simplify(y1*y2*y3-y1*y1-y3*y3*y0));
ygnum=ygN/ygD;
yg=simplify(ygnum/ybnum)   % yg - ygamma

%% 
disp('Linear stability with coupling but without diffusion')
y30=subs(y3,k,0)
y20=subs(y2,k,0)
y10=subs(y1,k,0)
y00=subs(y0,k,0)
yb0=subs(yb,k,0)
yg0=subs(yg,k,0)

%% Dispersion relation: 
clear Mr Mi
nn=100;
k0=linspace(0,kplus,nn)';
eta0=linspace(0,eta_end,nn)';

Mr=dispersion(CPy,k,eta,k0,eta0);
f1=figure; 
surf(k0,eta0,Mr,'linestyle','none');
xlabel('k')
ylabel('\eta')
title('Max(Re(\lambda))')
clim([-0.1,0.1])
view(2)
shading interp
ax = gca; 
ax.FontSize = 16;

cmap_colorbar()

saveas(f1,'DCy(k,eta).fig');
saveas(f1,'DCy(k,eta).png');

%% Surface Graphs: y3, y2, y1, y0, ybeta-yb, ygamma-yg
clear RHv
y0k=subs(y0,h,k^2); % y3 % y2 % y1 % yb % yg
nn=500;
k0=linspace(0,kplus,nn)';
eta0=linspace(0,eta_end,nn)';

[RHv,k00,eta00]=ycoeff(y0k,k0,eta0);

f2=figure; 
surf(k00,eta00,RHv,'linestyle','none')
xlabel('k')
ylabel('\eta')
zlabel('y_{0}(k,\eta)')
title('y_{0}(k,\eta)')
clim([-1,1])
% view(2)
shading interp
ax = gca; 
ax.FontSize = 18;

cmap_colorbar()

saveas(f2,'y_0(k,eta).fig'); % y3 % y2 % y1 % yb % yg
saveas(f2,'y_0(k,eta).png');

%% discriminant
y0h=subs(y0,k,sqrt(h));

[cfsy0h,hs5]=coeffs(y0h,h);
[discy0,P,D]=disc4(cfsy0h);
e_disc=double(solve(discy0==0,eta));
for i=1:length(e_disc)
p(i)=isreal(e_disc(i));
end
er=e_disc(p);
p1=find(e_disc(p)>0);
erp=er(p1);
disp('Roots of the discriminant are:')
disp(erp)
erp=[0 erp' 10];
erp=round(erp,2);
erpP=solve(P==0,eta);
erpD=solve(D==0,eta);
% figure; fplot(discy0,[0,10])
% yline(0,'b-','HandleVisibility', 'off');
% title('Discriminant of y0(eta)')

%% Dispersion relation: Scaled
clear Mr Mi
nn=100;
k0=linspace(0,kplus,nn)';
e1=linspace(0,0.17,0.1*nn);
e2=linspace(0.17,0.43,0.1*nn+1);
e3=linspace(0.43,1.41,0.2*nn+1);
e4=linspace(1.41,3.77,0.2*nn+1);
e5=linspace(3.77,10,0.4*nn+1);
eta0s=[e1 e2(2:end) e3(2:end) e4(2:end) e5(2:end)]'; % 

Mr=dispersion(CPy,k,eta,k0,eta0s);
axx=[0,0.17,0.43,1.41,3.77,10]; 
f3=figure; 
surf(k0,categorical(eta0s),Mr,'linestyle','none');
xlabel('k')
ylabel('\eta')
title('Max(Re(\lambda))')
clim([-0.1,0.1])
view(2)
shading interp
set(gca,'ytick',categorical(axx))
ax = gca; 
ax.FontSize = 16;

yline(categorical(0.17),'r--','Linewidth',2)
yline(categorical(0.43),'r--','Linewidth',2)
yline(categorical(1.41),'r--','Linewidth',2)
yline(categorical(3.77),'r--','Linewidth',2)

cmap_colorbar()

saveas(f3,'DCy(k,eta)_shrunk.fig');
saveas(f3,'DCy(k,eta)_shrunk.png');

%% y_0(k,eta): Scaled
clear RHv
y0k=subs(y0,h,k^2);
nn=500;
k0=linspace(0,kplus,nn)';
e1=linspace(0,0.17,0.1*nn);
e2=linspace(0.17,0.43,0.1*nn+1);
e3=linspace(0.43,1.41,0.2*nn+1);
e4=linspace(1.41,3.77,0.2*nn+1);
e5=linspace(3.77,10,0.4*nn+1);
eta0s=[e1 e2(2:end) e3(2:end) e4(2:end) e5(2:end)]'; % 

[RHv,k00,eta00]=ycoeff(y0k,k0,eta0s);

f4=figure; 
surf(k00,categorical(eta00),RHv,'linestyle','none')
xlabel('k')
ylabel('\eta')
zlabel('y_{0}(k,\eta)')
title('y_{0}(k,\eta)')
clim([-1,1])
view(2)
shading interp
axx=[0,0.17,0.43,1.41,3.77,10];
set(gca,'ytick',categorical(axx))
ax = gca; 
ax.FontSize = 18;

yline(categorical(0.17),'r--','Linewidth',2)
yline(categorical(0.43),'r--','Linewidth',2)
yline(categorical(1.41),'r--','Linewidth',2)
yline(categorical(3.77),'r--','Linewidth',2)

cmap_colorbar()

saveas(f4,'y_0(k,eta)_shrunk.fig');
saveas(f4,'y_0(k,eta)_shrunk.png');

%% Dispersion curves at different \eta values
clear Mr Mi
k0=linspace(0,kplus,100)';
for j=1:length(erp)
CPe=subs(CPy,eta,erp(j));
for i=1:length(k0)
    CP1=subs(CPe,k,k0(i));
    EG1=solve(CP1);
    rp=real(EG1);
    img=imag(EG1);
    Mr(i)=max(rp);
    % Mi(i)=max(img);
end
Mr=double(Mr);
f5=figure;
hold on
plot(k0,Mr,'Linewidth',3);
yline(0,'r--','HandleVisibility', 'off','Linewidth',3);
xlabel('k','FontSize',16);
ylabel('max(Re(\lambda))','FontSize',16);
title(['\eta=',num2str(erp(j))],'FontSize',16)
ylim([-1.5,0.5])
xlim([0,1])
ax = gca; 
ax.FontSize = 16;
saveas(f5,['eta=',num2str(erp(j)),'.fig']);
saveas(f5,['eta=',num2str(erp(j)),'.png']);
end

%% Reduced polynomials
y0_4dh=sum(cfsy0h(2:end).*hs5(2:end));
% del4=disc3(cfsy0(2:end));
% solve(del4==0)

y0_3dh=sum(cfsy0h(1,3:end).*hs5(1,3:end));
% del3=disc2(cfsy0(3:end));
% solve(del3==0)

y0_43dh=sum(cfsy0h(3:end).*hs5(3:end));
% del43=disc2(cfsy0(3:end));
% solve(del43==0)

%% y_0(k,eta): reduced y0_4d, y0_3d, y0_43d
clear RHv
y0k=subs(y0_4dh,h,k^2); % y0_3dh % y0_43dh
nn=500;
k0=linspace(0,kplus,nn)';
e1=linspace(0,0.17,0.1*nn);
e2=linspace(0.17,0.43,0.1*nn+1);
e3=linspace(0.43,1.41,0.2*nn+1);
e4=linspace(1.41,3.77,0.2*nn+1);
e5=linspace(3.77,10,0.4*nn+1);
eta0s=[e1 e2(2:end) e3(2:end) e4(2:end) e5(2:end)]'; % 

[RHv,k00,eta00]=ycoeff(y0k,k0,eta0s);

f6=figure; 
surf(k00,categorical(eta00),RHv,'linestyle','none')
xlabel('k')
ylabel('\eta')
% zlabel('y_{r_0}(k,\eta)')
% title('y_{r_0}(k,\eta)')
clim([-1,1])
view(2)
shading interp
axx=[0,0.17,0.43,1.41,3.77,10]; %
set(gca,'ytick',categorical(axx))
ax = gca; 
ax.FontSize = 18;

yline(categorical(0.17),'r--','Linewidth',2)
yline(categorical(0.43),'r--','Linewidth',2)
yline(categorical(1.41),'r--','Linewidth',2)
yline(categorical(3.77),'r--','Linewidth',2)

cmap_colorbar()

saveas(f6,'y_0_4d(k,eta).fig'); % y0_3d % y0_43d
saveas(f6,'y_0_4d(k,eta).png'); % y0_3d % y0_43d

%% Comparison Graph
nn=250;
k0=linspace(0,kplus,nn)';

sss=y0_red(y0h,h,k,eta,k0,10);
f7=figure(100);
hold on
plot(k0,sss,'o','LineWidth',3)

sss=y0_red(y0_4dh,h,k,eta,k0,10);
figure(100);
hold on
plot(k0,sss,'--','LineWidth',3)

sss=y0_red(y0_3dh,h,k,eta,k0,10);
figure(100);
hold on
plot(k0,sss,'o','LineWidth',3)

sss=y0_red(y0_43dh,h,k,eta,k0,10);
figure(100);
hold on
plot(k0,sss,'--','LineWidth',3)

ylim([-1000,1000])
xlabel('k')
% ylabel('y_{r_0}')
title('\eta=10')
xlim([0,1])
yline(0,'k--','Linewidth',2,'HandleVisibility', 'off')
legend('y_0','Removed 4^{th} degree', ...
    'Removed 3^{rd} degree','Removed both 4^{th} and 3^{rd} degree')
ax = gca; 
ax.FontSize = 16;
saveas(f7,'Compare_y0.fig');
saveas(f7,'Compare_y0.png');

%% Dispersion relation: Reduced
CPyh=subs(CPy,k,sqrt(h));
[cfsCPyh,yCPyyh]=coeffs(CPyh,y);
cfsCPyh(end)=y0_43dh;
CPred=sum(cfsCPyh.*yCPyyh);
CPredk=subs(CPred,h,k^2);
disp('Reduced Dispersion relation')
disp(CPredk)

% clear Mr
% nn=100;
% k0=linspace(0,kplus,nn)';
% eta0=linspace(0,eta_end,nn)';
% 
% Mr=dispersion(CPredk,k,eta,k0,eta0);
% 
% f2=figure; 
% surf(k0,eta0,Mr,'linestyle','none');
% xlabel('k')
% ylabel('\eta')
% title('Reduced Dispersion relation')
% clim([-0.1,0.1])
% view(2)
% shading interp
% ax = gca; 
% ax.FontSize = 16;
% cmap_colorbar()

%%
fprintf('\nDone!\n');
toc
diary OFF
