%% norme (k=1, sinsin, quadrato unitario, uniforme)
clear all; close all; clc;
h=1./[8; 16; 32; 64; 100];
l2=[0.00950273; 0.00267849; 0.000694077; 0.000175147; 7.18731e-5];
h1=[0.0836087; 0.0239605; 0.00618177; 0.00155723; 0.000638804];
infin=[0.0200595; 0.00548022; 0.00139685; 0.000350852; 0.000143841];

%% norme (k=1, sinsin, quadrato unitario, NONuniforme)
clear all; close all; clc;
h=[0.354802; 0.19089; 0.092233; 0.050394; 0.0243531];
l2=[0.912422; 0.0229184; 0.00437779; 0.00112855; 0.000281087];
h1=[0.398415; 0.194524; 0.0721927; 0.0363718; 0.0176535];
infin=[0.118301; 0.0534065; 0.0135009; 0.00341337; 0.00117141];

%% norme (k=2, sinsin, quadrato unitario, uniforme)
clear all; close all; clc;
h=1./[4; 8; 16; 32];
l2=[0.00913455; 0.000726485; 6.17313e-5; 4.14001e-6];
h1=[0.0602427; 0.00919023; 0.00145336; 0.000193643];
infin=[0.0156314; 0.000813505; 6.46555e-5; 4.28063e-6];

%% norme (k=4, sinsin, quadrato unitario, uniforme)
clear all; close all; clc;
h=1./[4; 8; 16; 32];
l2=[0.00502933; 9.36046e-5; 1.99072e-6; 3.1243e-8];
h1=[0.0229819; 0.000910436; 4.07612e-5; 1.2736e-6];
infin=[0.00269752; 5.35323e-5; 1.24377e-6; 2.7481e-8];


%% norme (elliptic, k=1, sinsin, quadrato unitario, uniforme)
clear all; close all; clc;
h=1./[8; 16; 32; 64; 100];
l2=[0.0125821; 0.00338597; 0.000867008; 0.000218131; 8.94594e-5];
h1=[0.102817; 0.0281694; 0.00719932; 0.00180968; 0.000742041];
infin=[0.0345832; 0.00901645; 0.00227893; 0.000573065; 0.000234797];

%% norme (elliptic, k=1, sinsin, quadrato unitario, NONuniforme)
clear all; close all; clc;
h=[0.354802; 0.19089; 0.092233; 0.050394; 0.0243531];
l2=[0.0989859; 0.026812; 0.00474022; 0.00118312; 0.000295651];
h1=[0.459277; 0.224952; 0.0746551; 0.0369677; 0.0176951];
infin=[0.132657; 0.0650969; 0.0155406; 0.00415896; 0.00125053];

%% norme (elliptic, k=4, sinsin, quadrato unitario, uniforme)
clear all; close all; clc;
h=1./[4; 8; 16; 32];
l2=[0.00419042; 8.91735e-5; 2.0202e-6; 3.1566e-8];
h1=[0.018986;  0.000848244; 4.11333e-5; 1.2854e-6];
infin=[0.00234457; 6.88552e-5; 1.68693e-6; 5.17562e-8];

%% norme (elliptic, k=2, sinsin, quadrato unitario, uniforme)
clear all; close all; clc;
h=1./[4; 8; 16; 32];
l2=[0.0256206; 0.00312888; 0.000234247; 1.5342e-5];
h1=[0.12818; 0.0274897; 0.00413025; 0.000544367];
infin=[0.0334666; 0.00348365; 0.000279229; 1.86031e-5];


%% norme (elliptic, k=2, sinsin, quadrato unitario, NONuniforme)
clear all; close all; clc;
h=[0.354802; 0.19089; 0.092233; 0.050394];
l2=[0.0394522; 0.00411537; 0.000621004; 6.78268e-5];
h1=[0.173823; 0.036949; 0.0104979; 0.00232494];
infin=[0.0403815; 0.0061549; 0.000945599; 0.000113993];

%%
figure;
loglog(h,l2,'bo-','linewidth',2);
hold on; grid on; box on;
esponente=3; esponente_super=4;
costante=l2(1)/h(1)^esponente_super-1; costante2=l2(1)/h(1)^esponente+2;
loglog(h,costante2*h.^esponente,'r','linewidth',2);
loglog(h,costante*h.^esponente_super,'g','linewidth',2);

xlabel({'$$h$$'},'interpreter','latex','fontsize',12);
ylabel({'$$||u-u_h||_{L^2(\Omega)}$$'},'interpreter','latex','fontsize',12);
%legend({'Numerical','$$h^3$$'},'interpreter','latex','fontsize',12);
legend({'Numerical','$$h^3$$','$$h^4$$'},'interpreter','latex','fontsize',12);
title({'Discretization error in $$L^2$$-norm'},'interpreter','latex','fontsize',12);
stima=log(l2(2:end)./l2(1:end-1))./log(h(2:end)./h(1:end-1))
stima_avg=mean(stima)
%%
figure;
loglog(h,h1,'bo-','linewidth',2);
hold on; grid on; box on;
esponente=2; esponente_super=3;
costante=h1(1)/h(1)^esponente_super-5; costante2=h1(1)/h(1)^esponente+1;
loglog(h,costante2*h.^esponente,'r','linewidth',2);
%loglog(h,costante*h.^esponente_super,'g','linewidth',2);

xlabel({'$$h$$'},'interpreter','latex','fontsize',12);
ylabel({'$$||u-u_h||_{H^1_0(\Omega)}$$'},'interpreter','latex','fontsize',12);
%legend({'Numerical','$$h$$'},'interpreter','latex','fontsize',12);
legend({'Numerical','$$h^2$$','$$h^3$$'},'interpreter','latex','fontsize',12);
title({'Discretization error in $$H_0^1$$-norm'},'interpreter','latex','fontsize',12);
stima=log(h1(2:end)./h1(1:end-1))./log(h(2:end)./h(1:end-1))
stima_avg=mean(stima)

%%
figure;
loglog(h,infin,'bo-','linewidth',2);
hold on; grid on; box on;
esponente=2; esponente_super=3;
costante=infin(1)/h(1)^esponente_super+2; costante2=infin(1)/h(1)^esponente+1;
loglog(h,costante2*h.^esponente,'r','linewidth',2);
loglog(h,costante*h.^esponente_super,'g','linewidth',2);

xlabel({'$$h$$'},'interpreter','latex','fontsize',12);
ylabel({'$$||u-u_h||_{L^\infty(\Omega)}$$'},'interpreter','latex','fontsize',12);
%legend({'Numerical','$$h^2$$'},'interpreter','latex','fontsize',12);
legend({'Numerical','$$h^2$$','$$h^3$$'},'interpreter','latex','fontsize',12);
title({'Discretization error in $$L^\infty$$-norm'},'interpreter','latex','fontsize',12);
stima=log(infin(2:end)./infin(1:end-1))./log(h(2:end)./h(1:end-1))
stima_avg=mean(stima)