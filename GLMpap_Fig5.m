% Continuum Model for Neurite Outgrowth
% Graham, Lauchlan & McLean Figure 5
% Steady-state length for variations in D, a and g
% (calculated via steady-state analysis)
% Version 2.0 (BPG & DRM 10-8-05)

% Parameters

% simulation (not used for analysis)
simp.dt = 0.01;                % time step
simp.tmax = 2000;              % simulation time
simp.datat = 100;              % data collection time step
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
k1 = 0.5;
k2 = 0.00001;
modp.e0 = modp.g*modp.sg/(k1*modp.c0*modp.rg*modp.a);  % soma flux-source rate
theta = 0;                     % fractional autoregulation
modp.er = theta*modp.e0;       % soma tubulin autoregulation
modp.rdt = 0;                  % autoregulation time delay
modp.el = k2*modp.rg;          % growth cone flux-sink rate
modp.zl = k2*modp.sg;          % growth cone flux-source rate

% plot parameters
tfs = 12;   % title font size

% Range for diffusion (D), active transport (a) and decay (g)
D = [30000 25000 20000 15000 10000 5000];
a = [100 80 60 40 20];
g = [0.002 0.0025 0.003 0.0035 0.005 0.01];

% Large growth
k1 = 0.5;
modp.e0 = modp.g*modp.sg/(k1*modp.c0*modp.rg*modp.a);  % soma flux-source rate
for (i=1:length(D))
  modp.D = D(i)
  % calculated parameters
  [calcp] = CMNG_calcparams(simp, modp);
  % get analytical steady-state values (note: not valid for ah=1)
  [CinflD, linflD(i)] = CMNG_lCanal(simp, modp, calcp, 0);
  linflD(i) = linflD(i)*(modp.D/(modp.rg*modp.c0));
  CinflD = CinflD*modp.c0;
  ahl = calcp.gamma*calcp.beta/(calcp.phi*calcp.alpha)
end
modp.D = 30000;
for (i=1:length(a))
  modp.a = a(i)
  % calculated parameters
  [calcp] = CMNG_calcparams(simp, modp);
  % get analytical steady-state values (note: not valid for ah=1)
  [Cinfla, linfla(i)] = CMNG_lCanal(simp, modp, calcp, 0);
  linfla(i) = linfla(i)*(modp.D/(modp.rg*modp.c0))
  Cinfla = Cinfla*modp.c0;
  ahl = calcp.gamma*calcp.beta/(calcp.phi*calcp.alpha)
end
modp.a = 100;
for (i=1:length(g))
  modp.g = g(i)
  % calculated parameters
  [calcp] = CMNG_calcparams(simp, modp);
  % get analytical steady-state values (note: not valid for ah=1)
  [Cinflg, linflg(i)] = CMNG_lCanal(simp, modp, calcp, 0);
  linflg(i) = linflg(i)*(modp.D/(modp.rg*modp.c0))
  Cinflg = Cinflg*modp.c0;
  ahl = calcp.gamma*calcp.beta/(calcp.phi*calcp.alpha)
end

% Moderate growth
k1 = 1.05;  % analysis not valid for k1=1
modp.e0 = modp.g*modp.sg/(k1*modp.c0*modp.rg*modp.a);  % soma flux-source rate
for (i=1:length(D))
  modp.D = D(i)
  % calculated parameters
  [calcp] = CMNG_calcparams(simp, modp);
  % get analytical steady-state values (note: not valid for ah=1)
  [CinfmD, linfmD(i)] = CMNG_lCanal(simp, modp, calcp, 0);
  linfmD(i) = linfmD(i)*(modp.D/(modp.rg*modp.c0));
  CinfmD = CinfmD*modp.c0;
  ahl = calcp.gamma*calcp.beta/(calcp.phi*calcp.alpha)
end
modp.D = 30000;
for (i=1:length(a))
  modp.a = a(i)
  % calculated parameters
  [calcp] = CMNG_calcparams(simp, modp);
  % get analytical steady-state values (note: not valid for ah=1)
  [Cinfma, linfma(i)] = CMNG_lCanal(simp, modp, calcp, 0);
  linfma(i) = linfma(i)*(modp.D/(modp.rg*modp.c0))
  Cinfma = Cinfma*modp.c0;
  ahl = calcp.gamma*calcp.beta/(calcp.phi*calcp.alpha)
end
modp.a = 100;
for (i=1:length(g))
  modp.g = g(i)
  % calculated parameters
  [calcp] = CMNG_calcparams(simp, modp);
  % get analytical steady-state values (note: not valid for ah=1)
  [Cinfmg, linfmg(i)] = CMNG_lCanal(simp, modp, calcp, 0);
  linfmg(i) = linfmg(i)*(modp.D/(modp.rg*modp.c0))
  Cinfmg = Cinfmg*modp.c0;
  ahl = calcp.gamma*calcp.beta/(calcp.phi*calcp.alpha)
end

% Small growth
k1 = 10;
modp.e0 = modp.g*modp.sg/(k1*modp.c0*modp.rg*modp.a);  % soma flux-source rate
for (i=1:length(D))
  modp.D = D(i)
  % calculated parameters
  [calcp] = CMNG_calcparams(simp, modp);
  % get analytical steady-state values (note: not valid for ah=1)
  [CinfsD, linfsD(i)] = CMNG_lCanal(simp, modp, calcp, 0);
  linfsD(i) = linfsD(i)*(modp.D/(modp.rg*modp.c0));
  CinfsD = CinfsD*modp.c0;
  ahl = calcp.gamma*calcp.beta/(calcp.phi*calcp.alpha)
end
modp.D = 30000;
for (i=1:length(a))
  modp.a = a(i)
  % calculated parameters
  [calcp] = CMNG_calcparams(simp, modp);
  % get analytical steady-state values (note: not valid for ah=1)
  [Cinfsa, linfsa(i)] = CMNG_lCanal(simp, modp, calcp, 0);
  linfsa(i) = linfsa(i)*(modp.D/(modp.rg*modp.c0))
  Cinfsa = Cinfsa*modp.c0;
  ahl = calcp.gamma*calcp.beta/(calcp.phi*calcp.alpha)
end
modp.a = 100;
for (i=1:length(g))
  modp.g = g(i)
  % calculated parameters
  [calcp] = CMNG_calcparams(simp, modp);
  % get analytical steady-state values (note: not valid for ah=1)
  [Cinfsg, linfsg(i)] = CMNG_lCanal(simp, modp, calcp, 0);
  linfsg(i) = linfsg(i)*(modp.D/(modp.rg*modp.c0))
  Cinfsg = Cinfsg*modp.c0;
  ahl = calcp.gamma*calcp.beta/(calcp.phi*calcp.alpha)
end

% Plot sensitivity against base case (D=30000, a=100, g=0.002)
dD = D./D(1);
da = a./a(1);
dg = g(1)./g;

subplot(3,3,1);
plot(dD, linflD./linflD(1), 'k*-');
hold on;
plot(dD, linfmD./linfmD(1), 'ko-');
plot(dD, linfsD./linfsD(1), 'ks-');
title('D','FontSize',tfs);
xlabel('Rate (%)');
ylabel('Length (%)');
%legend('large','mod','small',4);

subplot(3,3,2);
plot(da, linfla./linfla(1), 'k*-');
hold on;
plot(da, linfma./linfma(1), 'ko-');
plot(da, linfsa./linfsa(1), 'ks-');
title('a','FontSize',tfs);
xlabel('Rate (%)');

subplot(3,3,3);
plot(dg, linflg./linflg(1), 'k*-');
hold on;
plot(dg, linfmg./linfmg(1), 'ko-');
plot(dg, linfsg./linfsg(1), 'ks-');
title('g','FontSize',tfs);
xlabel('Rate (%)');
%legend('large','mod','small', 2);
