% Continuum Model for Neurite Outgrowth
% Graham, Lauchlan & McLean Figure 4
% Variations in D, a and g for large, medium and small growth regimes
%  - length profiles
% Version 1.0 (BPG & DRM 5-2-05)

% Parameters

% simulation
simp.dt = 0.01;                % time step
simp.tmax = 5000;              % simulation time
simp.datat = 1000;             % data collection time step
simp.N = 100;                  % number of spatial points
simp.kmax = 10000;             % maximum corrector steps
simp.mc = 0.0001;              % tolerance on C;
simp.ml = 0.0001;              % tolerance on l;

% user-defined
modp.c0 = 10;                  % concentration scale
modp.l0 = 0.01;                % initial (min) length;
modp.D = 30000;                % diffusion constant
modp.a = 100;                  % active transport rate
modp.g = 0.002;                % decay rate
modp.rg = 10;                  % growth rate constant
modp.sg = 100;                 % growth rate set point (threshold)
k1 = 0.5;                      % alpha_twid_h value
k2 = 0.00001;                  % assembly to concentration scale
modp.e0 = modp.g*modp.sg/(k1*modp.c0*modp.rg*modp.a);  % soma flux-source rate
theta = 0;                     % fractional autoregulation
modp.er = theta*modp.e0;       % soma tubulin autoregulation
modp.rdt = 0;                  % autoregulation time delay
modp.el = k2*modp.rg;          % growth cone flux-sink rate
modp.zl = k2*modp.sg;          % growth cone flux-source rate

% plot parameters
tfs = 12;   % title font size


% Run simulations

% Large growth regime
% Run 1: D=30000, a=100, g=0.002
% calculated parameters
[calcp] = CMNG_calcparams(simp, modp);
% run model for jmax time steps, linear ICs, no retraction
[Cl, C0l, CNl, ll] = CMNG_run(simp, modp, calcp, -1, modp);
[t, Cl, C0l, CNl, ll] = CMNG_dimen(simp, modp, Cl, C0l, CNl, ll);  % dimensionalise
Cal = [C0l Cl CNl];
% Run D2: D=20000
modp.D=20000;
% calculated parameters
[calcp] = CMNG_calcparams(simp, modp);
% run model for jmax time steps, linear ICs, no retraction
[ClD2, C0lD2, CNlD2, llD2] = CMNG_run(simp, modp, calcp, -1, modp);
[tD2, ClD2, C0lD2, CNlD2, llD2] = CMNG_dimen(simp, modp, ClD2, C0lD2, CNlD2, llD2);  % dimensionalise
% Run D3: D=10000
modp.D=10000;
% calculated parameters
[calcp] = CMNG_calcparams(simp, modp);
% run model for jmax time steps, linear ICs, no retraction
[ClD3, C0lD3, CNlD3, llD3] = CMNG_run(simp, modp, calcp, -1, modp);
[tD3, ClD3, C0lD3, CNlD3, llD3] = CMNG_dimen(simp, modp, ClD3, C0lD3, CNlD3, llD3);  % dimensionalise
% Run a2: a=80
modp.D=30000;
modp.a=80;
% calculated parameters
[calcp] = CMNG_calcparams(simp, modp);
% run model for jmax time steps, linear ICs, no retraction
[Cla2, C0la2, CNla2, lla2] = CMNG_run(simp, modp, calcp, -1, modp);
[t, Cla2, C0la2, CNla2, lla2] = CMNG_dimen(simp, modp, Cla2, C0la2, CNla2, lla2);  % dimensionalise
% Run a3: a=60
modp.a=60;
% calculated parameters
[calcp] = CMNG_calcparams(simp, modp);
% run model for jmax time steps, linear ICs, no retraction
[Cla3, C0la3, CNla3, lla3] = CMNG_run(simp, modp, calcp, -1, modp);
[t, Cla3, C0la3, CNla3, lla3] = CMNG_dimen(simp, modp, Cla3, C0la3, CNla3, lla3);  % dimensionalise
% Run g2: g=0.0025
modp.a=100;
modp.g=0.002*100/80;
% calculated parameters
[calcp] = CMNG_calcparams(simp, modp);
% run model for jmax time steps, linear ICs, no retraction
[Clg2, C0lg2, CNlg2, llg2] = CMNG_run(simp, modp, calcp, -1, modp);
[t, Clg2, C0lg2, CNlg2, llg2] = CMNG_dimen(simp, modp, Clg2, C0lg2, CNlg2, llg2);  % dimensionalise
% Run g3: g=0.0033
modp.g=0.002*100/60;
% calculated parameters
[calcp] = CMNG_calcparams(simp, modp);
% run model for jmax time steps, linear ICs, no retraction
[Clg3, C0lg3, CNlg3, llg3] = CMNG_run(simp, modp, calcp, -1, modp);
[t, Clg3, C0lg3, CNlg3, llg3] = CMNG_dimen(simp, modp, Clg3, C0lg3, CNlg3, llg3);  % dimensionalise

% Plot results
subplot(3,3,1);
plot(t,ll,'k-');
hold on;
plot(tD2,llD2,'k--');
plot(tD3,llD3,'k-.');
title('Large','FontSize',tfs);
%xlabel('Time');
ylabel('Length (\mum)');
%legend('D=30000','D=15000','D=6000');
subplot(3,3,4);
plot(t,ll,'k-');
hold on;
plot(t,lla2,'k--');
plot(t,lla3,'k-.');
%title('Length');
%xlabel('Time');
ylabel('Length (\mum)');
%legend('a=100','a=80','a=60');
subplot(3,3,7);
plot(t,ll,'k-');
hold on;
plot(t,llg2,'k--');
plot(t,llg3,'k-.');
%title('Length');
xlabel('Time (hrs)');
ylabel('Length (\mum)');
%legend('g=0.002','g=0.0025','g=0.0033');


% Moderate growth regime
modp.g = 0.002;
k1 = 1;
modp.e0 = modp.g*modp.sg/(k1*modp.c0*modp.rg*modp.a);  % soma flux-source rate
% Run 1: D=30000, a=100, g=0.002
% calculated parameters
[calcp] = CMNG_calcparams(simp, modp);
% run model for jmax time steps, linear ICs, no retraction
[Cm, C0m, CNm, lm] = CMNG_run(simp, modp, calcp, -1, modp);
[t, Cm, C0m, CNm, lm] = CMNG_dimen(simp, modp, Cm, C0m, CNm, lm);  % dimensionalise
% Run D2: D=20000
modp.D=20000;
% calculated parameters
[calcp] = CMNG_calcparams(simp, modp);
% run model for jmax time steps, linear ICs, no retraction
[CmD2, C0mD2, CNmD2, lmD2] = CMNG_run(simp, modp, calcp, -1, modp);
[tD2, CmD2, C0mD2, CNmD2, lmD2] = CMNG_dimen(simp, modp, CmD2, C0mD2, CNmD2, lmD2);  % dimensionalise
% Run D3: D=10000
modp.D=10000;
% calculated parameters
[calcp] = CMNG_calcparams(simp, modp);
% run model for jmax time steps, linear ICs, no retraction
[CmD3, C0mD3, CNmD3, lmD3] = CMNG_run(simp, modp, calcp, -1, modp);
[tD3, CmD3, C0mD3, CNmD3, lmD3] = CMNG_dimen(simp, modp, CmD3, C0mD3, CNmD3, lmD3);  % dimensionalise
% Run a2: a=80
modp.D=30000;
modp.a=80;
% calculated parameters
[calcp] = CMNG_calcparams(simp, modp);
% run model for jmax time steps, linear ICs, no retraction
[Cma2, C0ma2, CNma2, lma2] = CMNG_run(simp, modp, calcp, -1, modp);
[t, Cma2, C0ma2, CNma2, lma2] = CMNG_dimen(simp, modp, Cma2, C0ma2, CNma2, lma2);  % dimensionalise
% Run a3: a=60
modp.a=60;
% calculated parameters
[calcp] = CMNG_calcparams(simp, modp);
% run model for jmax time steps, linear ICs, no retraction
[Cma3, C0ma3, CNma3, lma3] = CMNG_run(simp, modp, calcp, -1, modp);
[t, Cma3, C0ma3, CNma3, lma3] = CMNG_dimen(simp, modp, Cma3, C0ma3, CNma3, lma3);  % dimensionalise
% Run g2: g=0.0025
modp.a=100;
modp.g=0.002*100/80;
% calculated parameters
[calcp] = CMNG_calcparams(simp, modp);
% run model for jmax time steps, linear ICs, no retraction
[Cmg2, C0mg2, CNmg2, lmg2] = CMNG_run(simp, modp, calcp, -1, modp);
[t, Cmg2, C0mg2, CNmg2, lmg2] = CMNG_dimen(simp, modp, Cmg2, C0mg2, CNmg2, lmg2);  % dimensionalise
% Run g3: g=0.0033
modp.g=0.002*100/60;
% calculated parameters
[calcp] = CMNG_calcparams(simp, modp);
% run model for jmax time steps, linear ICs, no retraction
[Cmg3, C0mg3, CNmg3, lmg3] = CMNG_run(simp, modp, calcp, -1, modp);
[t, Cmg3, C0mg3, CNmg3, lmg3] = CMNG_dimen(simp, modp, Cmg3, C0mg3, CNmg3, lmg3);  % dimensionalise

% Plot results
subplot(3,3,2);
plot(t,lm,'k-');
hold on;
plot(tD2,lmD2,'k--');
plot(tD3,lmD3,'k-.');
title('Moderate','FontSize',tfs);
%xlabel('Time');
%ylabel('Length');
legend('D=30','D=20','D=10');
subplot(3,3,5);
plot(t,lm,'k-');
hold on;
plot(t,lma2,'k--');
plot(t,lma3,'k-.');
%title('Length');
%xlabel('Time');
%ylabel('Length');
legend('a=100','a=80','a=60');
subplot(3,3,8);
plot(t,lm,'k-');
hold on;
plot(t,lmg2,'k--');
plot(t,lmg3,'k-.');
%title('Length');
xlabel('Time (hrs)');
%ylabel('Length');
legend('g=20','g=25','g=33');


% Small growth regime
simp.tmax = 200;                  % simulation time
simp.datat = 100;                 % data collection time step
modp.g = 0.002;
k1 = 10;
modp.e0 = modp.g*modp.sg/(k1*modp.c0*modp.rg*modp.a);  % soma flux-source rate
% Run 1: D=30000, a=100, g=0.002
% calculated parameters
[calcp] = CMNG_calcparams(simp, modp);
% run model for jmax time steps, linear ICs, no retraction
[Cs, C0s, CNs, ls] = CMNG_run(simp, modp, calcp, -1, modp);
[t, Cs, C0s, CNs, ls] = CMNG_dimen(simp, modp, Cs, C0s, CNs, ls);  % dimensionalise
% Run D2: D=20000
modp.D=20000;
% calculated parameters
[calcp] = CMNG_calcparams(simp, modp);
% run model for jmax time steps, linear ICs, no retraction
[CsD2, C0sD2, CNsD2, lsD2] = CMNG_run(simp, modp, calcp, -1, modp);
[tD2, CsD2, C0sD2, CNsD2, lsD2] = CMNG_dimen(simp, modp, CsD2, C0sD2, CNsD2, lsD2);  % dimensionalise
% Run D3: D=10000
modp.D=10000;
% calculated parameters
[calcp] = CMNG_calcparams(simp, modp);
% run model for jmax time steps, linear ICs, no retraction
[CsD3, C0sD3, CNsD3, lsD3] = CMNG_run(simp, modp, calcp, -1, modp);
[tD3, CsD3, C0sD3, CNsD3, lsD3] = CMNG_dimen(simp, modp, CsD3, C0sD3, CNsD3, lsD3);  % dimensionalise
% Run a2: a=80
modp.D=30000;
modp.a=80;
% calculated parameters
[calcp] = CMNG_calcparams(simp, modp);
% run model for jmax time steps, linear ICs, no retraction
[Csa2, C0sa2, CNsa2, lsa2] = CMNG_run(simp, modp, calcp, -1, modp);
[t, Csa2, C0sa2, CNsa2, lsa2] = CMNG_dimen(simp, modp, Csa2, C0sa2, CNsa2, lsa2);  % dimensionalise
% Run a3: a=60
modp.a=60;
% calculated parameters
[calcp] = CMNG_calcparams(simp, modp);
% run model for jmax time steps, linear ICs, no retraction
[Csa3, C0sa3, CNsa3, lsa3] = CMNG_run(simp, modp, calcp, -1, modp);
[t, Csa3, C0sa3, CNsa3, lsa3] = CMNG_dimen(simp, modp, Csa3, C0sa3, CNsa3, lsa3);  % dimensionalise
% Run g2: g=0.0025
modp.a=100;
modp.g=0.002*100/80;
% calculated parameters
[calcp] = CMNG_calcparams(simp, modp);
% run model for jmax time steps, linear ICs, no retraction
[Csg2, C0sg2, CNsg2, lsg2] = CMNG_run(simp, modp, calcp, -1, modp);
[t, Csg2, C0sg2, CNsg2, lsg2] = CMNG_dimen(simp, modp, Csg2, C0sg2, CNsg2, lsg2);  % dimensionalise
% Run g3: g=0.0033
modp.g=0.002*100/60;
% calculated parameters
[calcp] = CMNG_calcparams(simp, modp);
% run model for jmax time steps, linear ICs, no retraction
[Csg3, C0sg3, CNsg3, lsg3] = CMNG_run(simp, modp, calcp, -1, modp);
[t, Csg3, C0sg3, CNsg3, lsg3] = CMNG_dimen(simp, modp, Csg3, C0sg3, CNsg3, lsg3);  % dimensionalise

% Small growth
subplot(3,3,3);
plot(t,ls,'k-');
hold on;
plot(tD2,lsD2,'k--');
plot(tD3,lsD3,'k-.');
title('Small','FontSize',tfs);
%xlabel('Time');
%ylabel('Length');
%legend('D=30000','D=15000','D=6000');
subplot(3,3,6);
plot(t,ls,'k-');
hold on;
plot(t,lsa2,'k--');
plot(t,lsa3,'k-.');
%title('Length');
%xlabel('Time');
%ylabel('Length');
%legend('a=100','a=80','a=60');
subplot(3,3,9);
plot(t,ls,'k-');
hold on;
plot(t,lsg2,'k--');
plot(t,lsg3,'k-.');
%title('Length');
xlabel('Time (hrs)');
%ylabel('Length');
%legend('g=0.002','g=0.0025','g=0.0033');
