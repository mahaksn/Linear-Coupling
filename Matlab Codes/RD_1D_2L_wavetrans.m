clc; clear all; close all;
tic
% U-layer_1; L-layer_2;
% u,v - concentrations;

dim='1D';
layer='2L';
suffix='wavetrans';
pdfname='RD_1D_2L_wavetrans';

%% pre-setting
showanimation=1;
makegif=1;
GetPatterns=1;
drawperframe=200;
pattern=1; % Pattern fully developed at..

Schnakenberg=1;

case1=0;
case2=1; % non-dimensionalised

%% set the parameters which govern equation

T=5*10^3;
nx=200;
tol=1e-10;

tend=100;

if Schnakenberg
RD='Schnak';
f = @(u,v,a,b) (a - u + (u.^2).*v);
g = @(u,v,a,b) (b - (u.^2).*v);
if case1
    RD=[RD,'_','case1'];
    a=0.05; b=1.4;
    aU=a; bU=b; aL=a; bL=b;
    DuU=0.00025; DvU=0.0125;
    DuL=0.01; DvL=0.2;
    L=2;
    eta_end=10;
    kplus=60;
elseif case2
    RD=[RD,'_','case2'];
    a=0.05; b=1.4;
    aU=a; bU=b; aL=a; bL=b;
    DuU=1; DvU=50;
    DuL=40; DvL=800;
    L=130;
    eta_end=10;
    kplus=1;
end
end

%%
e1=linspace(0,0.17,nx/4);
e2=linspace(0.17,0.43,nx/4+1);
e3=linspace(0.43,3.77,nx/4+1);
e4=linspace(3.77,10,nx/4+1);
eta=[e1 e2(2:end) e3(2:end) e4(2:end)];

for e=1:length(eta)
eta_u = eta(e);
eta_v = eta_u; % fixed
pattern=1; % Pattern fully developed at..
tend=100;

%% saving
eta_text=['_eta=',num2str(eta(e))];
[filename,folder] = savingcode(RD,dim,suffix,layer);
prefix = strcat(filename,...
    eta_text);
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

%% domain setting
x=linspace(0,L,nx)';
dx=abs(x(2)-x(1));

%% time discretization

% T=500;%10^4;
dt=0.01;
nt=T/dt+1;
stopti=nt; % stop patterning at nt
nFrame=ceil((T/dt)/drawperframe);

%% FDM setup

o=ones(nx,1);
I=speye(nx,nx);
A=spdiags([o -2*o o],[-1 0 1],nx,nx);
A(1,2)=2; % for no-flux BC
A(nx,nx-1)=2;

% layer_1
CuU=DuU./(2*dx^2);
Lap_uU=CuU .* A;
CvU=DvU./(2*dx^2);
Lap_vU=CvU .* A;

% layer_2
CuL=DuL./(2*dx^2);
Lap_uL=CuL .* A;
CvL=DvL./(2*dx^2);
Lap_vL=CvL .* A;

% A-matrix
AuU = I - dt*Lap_uU;
AvU = I - dt*Lap_vU;
AuL = I - dt*Lap_uL;
AvL = I - dt*Lap_vL;

%% initial condition

pit=0; rng(pit); %change random seed
Perturbations = 0.01*(2*rand(nx,1)-1);

uU=uU00 + Perturbations;
vU=vU00 + Perturbations;
uL=uL00 + Perturbations;
vL=vL00 + Perturbations;

%% Numerical Simulation
if GetPatterns
%% Set up figure
giffile1 = [prefix,'_1','.gif'];
giffile2 = [prefix,'_2','.gif'];
if showanimation
    fig_pos = [100 100 1000 500];
    figU=figure('Position',fig_pos,'color','w');
    hold on
    % layer_1
    uUfig=plot(x,uU,'Linewidth',3);
    vUfig=plot(x,vU,'Linewidth',3);
    ylim([0,4]);
    xlabel('x');
    ylabel('u,v');
    legend('u_1','v_1')
    figtitle1=title('t=0');
    hold off
    
    figL=figure('Position',fig_pos,'color','w');
    hold on
    % layer_2
    uLfig=plot(x,uL,'Linewidth',3);
    vLfig=plot(x,vL,'Linewidth',3);
    ylim([0,4]);
    xlabel('x');
    ylabel('u,v');
    legend('u_2','v_2')
    figtitle2=title('t=0');
    hold off
end

%% simulation iteration

for ti=1:1:nt
    t=dt*(ti-1);
    if (mod(ti, drawperframe) == 1)
        if showanimation
            % layer_1
            uUfig.YData=uU;
            vUfig.YData=vU;
            % layer_2
            uLfig.YData=uL;
            vLfig.YData=vL;

            figtitle1.String=['t=',num2str(t,'%.1f')];
            figtitle2.String=['t=',num2str(t,'%.1f')];
            drawnow;
        end
        iFrame=(ti-1)/drawperframe+1;
        if makegif
            frame1 = getframe(figU);
            im1 = frame2im(frame1);
            [imind1,cm1] = rgb2ind(im1,256);
            if iFrame==1
                imwrite(imind1,cm1,giffile1,'gif', 'Loopcount',inf);
            else
                imwrite(imind1,cm1,giffile1,'gif','WriteMode', ...
                    'append','DelayTime',0);
            end
        end
        if makegif
            frame2 = getframe(figL);
            im2 = frame2im(frame2);
            [imind2,cm2] = rgb2ind(im2,256);
            if iFrame==1
                imwrite(imind2,cm2,giffile2,'gif', 'Loopcount',inf);
            else
                imwrite(imind2,cm2,giffile2,'gif','WriteMode', ...
                    'append','DelayTime',0);
            end
        end
    end
    
    cu = uL-uU;
    cv = vL-vU;

    % layer_1
    fU = f(uU,vU,aU,bU);
    gU = g(uU,vU,aU,bU);
    
    BuU = uU + dt*(Lap_uU*uU + fU + eta_u*cu);
    uUnew = AuU\BuU;
    BvU = vU + dt*(Lap_vU*vU + gU + eta_v*cv);
    vUnew = AvU\BvU;
    
    % layer_2
    fL = f(uL,vL,aL,bL);
    gL = g(uL,vL,aL,bL);
    
    BuL = uL + dt*(Lap_uL*uL + fL - eta_u*cu);
    uLnew = AuL\BuL;
    BvL = vL + dt*(Lap_vL*vL + gL - eta_v*cv);
    vLnew = AvL\BvL;
    
    PFD=[pattern*o,...
        (abs(uU-uUnew)<tol),...
        (abs(vU-vUnew)<tol),...
        (abs(uL-uLnew)<tol),...
        (abs(vL-vLnew)<tol)];

    if all(PFD)
        pattern=0;
        fprintf('\npattern fully developed at %.5f\n',t);
        tend=tend/dt+1;
        stopti=ti+tend;
    end
    
    if ti==stopti
        break
    end

    uU=uUnew; vU=vUnew;
    u1(e,:)=uUnew; v1(e,:)=vUnew;
    uL=uLnew; vL=vLnew;
    u2(e,:)=uLnew; v2(e,:)=vLnew;
end
end
%% saving final pattern

if showanimation
saveas(figU,[prefix,'_finalU.png']);
saveas(figU,[prefix,'_finalU.fig']);
saveas(figL,[prefix,'_finalL.png']);
saveas(figL,[prefix,'_finalL.fig']);
end
end

%%
axx=[0,0.17,0.43,3.77,10];

f1=figure; surf(x,categorical(eta),u1,'linestyle','none')
view(2)
xlabel('x')
ylabel('\eta')
title('u_1')
colorbar
set(gca,'ytick',categorical(axx))
ax = gca; 
ax.FontSize = 16;
saveas(f1,'u1_pattern_trans.fig');
saveas(f1,'u1_pattern_trans.png');

%%
f2=figure; surf(x,categorical(eta),v1,'linestyle','none')
view(2)
xlabel('x')
ylabel('\eta')
title('v_1')
colorbar
set(gca,'ytick',categorical(axx))
ax = gca; 
ax.FontSize = 16; 
saveas(f2,'v1_pattern_trans.fig');
saveas(f2,'v1_pattern_trans.png');

%%
f3=figure; surf(x,categorical(eta),u2,'linestyle','none')
view(2)
xlabel('x')
ylabel('\eta')
title('u_2')
colorbar
set(gca,'ytick',categorical(axx))
ax = gca; 
ax.FontSize = 16; 
saveas(f3,'u2_pattern_trans.fig');
saveas(f3,'u2_pattern_trans.png');

%%
f4=figure; surf(x,categorical(eta),v2,'linestyle','none')
view(2)
xlabel('x')
ylabel('\eta')
title('v_2')
colorbar
set(gca,'ytick',categorical(axx))
ax = gca; 
ax.FontSize = 16; 
saveas(f4,'v2_pattern_trans.fig');
saveas(f4,'v2_pattern_trans.png');

%%
fprintf('\nDone!\n');
toc
diary OFF
