% Continuum Model for Neurite Outgrowth with Autoregulation
% Graham, Lauchlan & McLean Figures 2 and 3
% Large, moderate & small growth regimes with dimensional parameters
% Version 1.0 (BPG & DRM 5-2-05)

% Parameters

% simulation
simp.dt = 0.01;                  % time step
simp.tmax = 5000;                % simulation time
simp.datat = 1000;               % data collection time step
simp.N = 100;                    % number of spatial points
simp.kmax = 10000;               % maximum corrector steps
simp.mc = 0.0001;                % tolerance on C;
simp.ml = 0.0001;                % tolerance on l;

% user-defined
modp.c0 = 10;                    % concentration scale
modp.l0 = 0.01;                  % initial (min) length;
modp.D = 30000;                  % diffusion constant
modp.a = 100;                    % active transport rate
modp.g = 0.002;                  % decay rate
modp.rg = 10;                    % growth rate constant
modp.sg = 100;                   % growth rate set point (threshold)
k1 = 0.5;                        % alpha_twid_h value
k2 = 0.00001;                    % assembly to concentration scale
modp.e0 = modp.g*modp.sg/(k1*modp.c0*modp.rg*modp.a);  % soma flux-source rate
theta = 0;                       % fractional autoregulation
modp.er = theta*modp.e0;         % soma tubulin autoregulation
modp.rdt = 0;                    % autoregulation time delay
modp.el = k2*modp.rg;            % growth cone flux-sink rate
modp.zl = k2*modp.sg;            % growth cone flux-source rate

% plot parameters
tfs = 12;   % title font size


doruns = 1;     % flag to run simulations

% Large growth
if (doruns == 1)
% Run 1: linear ICs
k1 = 0.5;
modp.e0 = modp.g*modp.sg/(k1*modp.c0*modp.rg*modp.a); 
% calculated parameters
[calcp] = CMNG_calcparams(simp, modp);
% run model for jmax time steps, linear ICs
[Cl1, C0l1, CNl1, ll1] = CMNG_run(simp, modp, calcp, -1, modp);
[t, Cl1, C0l1, CNl1, ll1] = CMNG_dimen(simp, modp, Cl1, C0l1, CNl1, ll1);  % dimensionalise
Cal1 = [C0l1 Cl1 CNl1];
% Run 2: initial small growth
modp.e0 = modp.e0/20;
% calculated parameters
[calcp] = CMNG_calcparams(simp, modp);
newp = modp;
newp.e0 = modp.e0*20;
% run model for jmax time steps, linear ICs, change to large growth
[Cl2, C0l2, CNl2, ll2] = CMNG_run(simp, modp, calcp, 100, newp);
[t, Cl2, C0l2, CNl2, ll2] = CMNG_dimen(simp, modp, Cl2, C0l2, CNl2, ll2);  % dimensionalise
Cal2 = [C0l2 Cl2 CNl2];
end

% Plot results (figure 1)
subplot(3,3,1);
plot(t,ll1,'k-');
hold on;
plot(t,ll2,'k-.');
title('Large','FontSize',tfs);
ylabel('Length (\mum)');
legend('linIC','grow',4);
subplot(3,3,4);
plot(t,C0l1,'k-');
hold on;
plot(t,C0l2,'k-.');
ylabel('Soma Concentration (\muM)');
subplot(3,3,7);
plot(t,CNl1,'k-');
hold on;
plot(t,CNl2,'k-.');
xlabel('Time (hours)');
ylabel('Terminal Concentration (\muM)');


% Moderate growth
if (doruns == 1)
% Run 1: linear ICs
k1 = 1;
modp.e0 = modp.g*modp.sg/(k1*modp.c0*modp.rg*modp.a); 
% calculated parameters
[calcp] = CMNG_calcparams(simp, modp);
% run model for jmax time steps, linear ICs
[Cm1, C0m1, CNm1, lm1] = CMNG_run(simp, modp, calcp, -1, modp);
[t, Cm1, C0m1, CNm1, lm1] = CMNG_dimen(simp, modp, Cm1, C0m1, CNm1, lm1);  % dimensionalise
Cam1 = [C0m1 Cm1 CNm1];
% Run 2: initial small growth
modp.e0 = modp.e0/10;
% calculated parameters
[calcp] = CMNG_calcparams(simp, modp);
newp = modp;
newp.e0 = modp.e0*10;
% run model for jmax time steps, linear ICs, change to moderate growth
[Cm2, C0m2, CNm2, lm2] = CMNG_run(simp, modp, calcp, 100, newp);
[t, Cm2, C0m2, CNm2, lm2] = CMNG_dimen(simp, modp, Cm2, C0m2, CNm2, lm2);  % dimensionalise
Cam2 = [C0m2 Cm2 CNm2];
end

% Plot results (figure 1)
subplot(3,3,2);
plot(t,lm1,'k-');
hold on;
plot(t,lm2,'k-.');
title('Moderate','FontSize',tfs);
subplot(3,3,5);
plot(t,C0m1,'k-');
hold on;
plot(t,C0m2,'k-.');
subplot(3,3,8);
plot(t,CNm1,'k-');
hold on;
plot(t,CNm2,'k-.');
xlabel('Time (hours)');


% Small growth
if (doruns == 1)
% Run 1: linear ICs
simp.tmax = 200;
simp.datat = 100;
k1 = 10;
modp.e0 = modp.g*modp.sg/(k1*modp.c0*modp.rg*modp.a); 
% calculated parameters
[calcp] = CMNG_calcparams(simp, modp);
% run model for jmax time steps, linear ICs
[Cs1, C0s1, CNs1, ls1] = CMNG_run(simp, modp, calcp, -1, modp);
[t, Cs1, C0s1, CNs1, ls1] = CMNG_dimen(simp, modp, Cs1, C0s1, CNs1, ls1);  % dimensionalise
Cas1 = [C0s1 Cs1 CNs1];
end

% Plot results (figure 1)
subplot(3,3,3);
plot(t,ls1,'k-');
hold on;
title('Small','FontSize',tfs);
subplot(3,3,6);
plot(t,C0s1,'k-');
hold on;
subplot(3,3,9);
plot(t,CNs1,'k-');
hold on;
xlabel('Time (hours)');


% Plot final concentration gradients (Figure 2)
% (assumes N does not change between runs)
figure;
space=0:1/simp.N:1;
subplot(3,3,1);
plot(space, Cal1(length(ll1),:), 'k-');
hold on;
plot(space, Cal2(length(ll2),:), 'k-.');
axis([0 1 10 20]);
title('Large','FontSize',tfs);
xlabel('Space');
ylabel('Concentration (\muM)');
legend('linIC','grow');
subplot(3,3,2);
plot(space, Cam1(length(lm1),:), 'k-');
hold on;
plot(space, Cam2(length(lm2),:), 'k-.');
axis([0 1 10 10.2]);
title('Moderate','FontSize',tfs);
xlabel('Space');
subplot(3,3,3);
plot(space, Cas1(length(ls1),:), 'k-');
hold on;
axis([0 1 10 10.00035]);
title('Small','FontSize',tfs);
xlabel('Space');

