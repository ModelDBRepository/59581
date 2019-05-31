% Continuum Model for Neurite Outgrowth with Autoregulation
% Example main m-file
% Simulation set up with dimensional parameters
% Version 1.0 (BPG & DRM 28-10-05)

% Parameters

% simulation
simp.dt = 0.01;                 % time step
simp.tmax = 100;                % simulation time
simp.datat = 100;               % data collection time step
simp.N = 100;                   % number of spatial points
simp.kmax = 10000;              % maximum corrector steps
simp.mc = 0.0001;               % tolerance on C;
simp.ml = 0.0001;               % tolerance on l;

% user-defined
modp.c0 = 10;                   % concentration scale
modp.l0 = 0.01;                 % initial (min) length;
modp.D = 30000;                 % diffusion constant
modp.a = 100;                   % active transport rate
modp.g = 0.002;                 % decay rate
modp.rg = 10;                   % growth rate constant
modp.sg = 100;                  % growth rate set point (threshold)
% Set alpha_th value 
% (<< 1 = large growth; around 1 = moderate; >> 1 = small)
alpha_th = 10;                  % alpha_twid_h value
modp.e0 = modp.g*modp.sg/(alpha_th*modp.c0*modp.rg*modp.a);  % soma flux-source rate
% Set tubulin production rate
% (alternative to setting alpha_th)
%modp.e0 = 0.000002;             % soma flux-source rate
%alpha_th = modp.g*modp.sg/(modp.e0*modp.c0*modp.rg*modp.a);  % alpha_twid_h value
theta = 0;                      % fractional autoregulation
modp.er = theta*modp.e0;        % soma tubulin autoregulation
modp.rdt = 0;                   % autoregulation time delay
growth_scale = 0.00001;
modp.el = growth_scale*modp.rg;           % growth cone flux-sink rate
modp.zl = growth_scale*modp.sg;           % growth cone flux-source rate

% Simulation with linear ICs
% calculated parameters
[calcp] = CMNG_calcparams(simp, modp);
% run model for jmax time steps
[C1, C01, CN1, l1] = CMNG_run(simp, modp, calcp, -1, modp);
[t, C1, C01, CN1, l1] = CMNG_dimen(simp, modp, C1, C01, CN1, l1);  % dimensionalise
Ca1 = [C01 C1 CN1];


% Plot results
plot3d = 1; % set to 1 if 3D concentration plot required

subplot(2,2,1);
plot(t,l1,'k-');
hold on;
title('Length');
xlabel('Time');
ylabel('Length');

subplot(2,2,2);
space=0:1/simp.N:1;
plot(space, Ca1(length(l1),:), 'k-');
hold on;
title('Concentration Gradient');
xlabel('Space');
ylabel('Concentration');

subplot(2,2,3);
plot(t,C01,'k-');
hold on;
title('Soma Concentration');
xlabel('Time');
ylabel('Concentration');

subplot(2,2,4);
plot(t,CN1,'k-');
hold on;
title('Terminal Concentration');
xlabel('Time');
ylabel('Concentration');


if (plot3d == 1)
figure(2);
surf(space,t,Ca1);
title('Concentration over Unit Space and Time');
xlabel('Space');
ylabel('Time');
zlabel('Concentration');
view(90-37.5,30);
end;

