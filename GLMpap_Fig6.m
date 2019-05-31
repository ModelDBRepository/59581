% Continuum Model for Neurite Outgrowth
% Graham, Lauchlan & McLean Figure 6
% Retraction from large to moderate growth regime
% (sg, rg or e0 change to cause retraction)
% Version 1.0 (BPG & DRM 7-2-05)

% Parameters

% simulation
simp.dt = 0.01;                % time step
simp.tmax = 10000;             % simulation time
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

doruns = 1;     % flag to run simulations
if (doruns == 1)

% Run 1: rg=5 half-way through
newp = modp;
newp.rg = 5;                    % growth rate set point (threshold)
newp.el = k2*newp.rg;           % growth cone flux-sink rate
[calcp] = CMNG_calcparams(simp, modp);  % calculated parameters
% run model for jmax time steps, linear ICs
[C1, C01, CN1, l1] = CMNG_run(simp, modp, calcp, 5000, newp);
%l1 = l1./2; % rescale for original rg
[t, C1, C01, CN1, l1] = CMNG_dimen(simp, modp, C1, C01, CN1, l1);  % dimensionalise
Ca1 = [C01 C1 CN1];
% get analytical steady-state values
[Cinfa1, linfa1] = CMNG_lCanal(simp, modp, calcp, 0);
linfa1 = linfa1*(modp.D/(modp.rg*modp.c0))
Cinfa1 = Cinfa1*modp.c0;
ah1 = calcp.gamma*calcp.beta/(calcp.phi*calcp.alpha)

% Run 2: sg=200 half-way through
modp.rg = 10;
modp.el = k2*modp.rg;            % growth cone flux-sink rate
newp = modp;
newp.sg = 200;                   % growth rate set point (threshold)
newp.zl = k2*newp.sg;            % growth cone flux-source rate
[calcp] = CMNG_calcparams(simp, modp);  % calculated parameters
% run model for jmax time steps, linear ICs
[C2, C02, CN2, l2] = CMNG_run(simp, modp, calcp, 5000, newp);
[t, C2, C02, CN2, l2] = CMNG_dimen(simp, modp, C2, C02, CN2, l2);  % dimensionalise
Ca2 = [C02 C2 CN2];
% get analytical steady-state values
[Cinfa2, linfa2] = CMNG_lCanal(simp, modp, calcp, 0);
linfa2 = linfa2*(modp.D/(modp.rg*modp.c0))
Cinfa2 = Cinfa2*modp.c0;
ah2 = calcp.gamma*calcp.beta/(calcp.phi*calcp.alpha)

% Run 3: e0 halved half-way through
modp.sg = 100;                    % growth rate set point (threshold)
modp.zl = k2*modp.sg;             % growth cone flux-source rate
newp = modp;
newp.e0=modp.e0/2;
[calcp] = CMNG_calcparams(simp, modp);  % calculated parameters
% run model for jmax time steps, linear ICs, no retraction
[C3, C03, CN3, l3] = CMNG_run(simp, modp, calcp, 5000, newp);
[t, C3, C03, CN3, l3] = CMNG_dimen(simp, modp, C3, C03, CN3, l3);  % dimensionalise
Ca3 = [C03 C3 CN3];
% get analytical steady-state values
[Cinfa3, linfa3] = CMNG_lCanal(simp, modp, calcp, 0);
linfa3 = linfa3*(modp.D/(modp.rg*modp.c0))
Cinfa3 = Cinfa3*modp.c0;
ah3 = calcp.gamma*calcp.beta/(calcp.phi*calcp.alpha)

end


% Plot results

subplot(3,3,1);
plot(t,l1,'k-');
hold on;
plot(t,l2,'k--');
plot(t,l3,'k-.');
title('Length','FontSize',tfs);
xlabel('Time');
ylabel('Length');
%legend('rg', 'sg','e0');

subplot(3,3,2);
plot(t,C01,'k-');
hold on;
plot(t,C02,'k--');
plot(t,C03,'k-.');
title('Soma Concentration','FontSize',tfs);
xlabel('Time');
ylabel('Concentration');
legend('rg', 'sg','e0');

subplot(3,3,3);
plot(t,CN1,'k-');
hold on;
plot(t,CN2,'k--');
plot(t,CN3,'k-.');
title('Terminal Concentration','FontSize',tfs);
xlabel('Time');
ylabel('Concentration');
