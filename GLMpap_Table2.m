% Continuum Model for Neurite Outgrowth
% Main m-file
% Comparison of steady-state lengths with Miller & Samuels, Table 1
% Version 1.0 (BPG & DRM 19-1-05)

% Parameters

% simulation
simp.dt = 0.01;                  % time step
simp.tmax = 2000;                  % simulation time
simp.datat = 100;                 % data collection time step
simp.N = 100;                     % number of spatial points
simp.kmax = 10000;               % maximum corrector steps
simp.mc = 0.0001;              % tolerance on C;
simp.ml = 0.0001;              % tolerance on l;

% user-defined
modp.c0 = 10;                     % concentration scale
modp.l0 = 0.01;                   % initial (min) length;
modp.D = 30000;                    % diffusion constant
modp.a = 100;                      % active transport rate
modp.g = 0.002;                    % decay rate
modp.rg = 10;                   % growth rate constant
modp.sg = 100;                   % growth rate set point (threshold)
k1 = 0.5;
k2 = 0.00001;
modp.e0 = modp.g*modp.sg/(k1*modp.c0*modp.rg*modp.a);                  % soma flux-source rate
modp.el = k2*modp.rg;                   % growth cone flux-sink rate
modp.zl = k2*modp.sg;                   % growth cone flux-source rate
%modp.g = 0;                    % decay rate


% Run 1: Garfish olfactory
modp.a = 2.38*1000/24;             % active transport rate (um/hr)
T05 = 75*24;                       % half-life (hours)
modp.g = 1/(T05/log(2));           % decay rate (per hour)
k1 = 0.34;                          % alpha_twid_h
modp.e0 = modp.g*modp.sg/(k1*modp.c0*modp.rg*modp.a);  % soma flux rate
% calculated parameters
[calcp] = CMNG_calcparams(simp, modp);
% run model for jmax time steps, linear ICs, no retraction
%[C1, C01, CN1, l1] = CMNG_rundyn(simp, modp, calcp, 0, 0);
%Ca1 = [C01 C1 CN1];
% get analytical steady-state values (note: not valid for ah=1)
[Cinfa1, linfa1] = CMNG_lCanal(simp, modp, calcp, 0);
linfa1 = linfa1*(modp.D/(modp.rg*modp.c0))/1000     % mm
Cinfa1 = Cinfa1*modp.c0;
ah1 = calcp.gamma*calcp.beta/(calcp.phi*calcp.alpha)
modp

% Run 2: Rat sciatic
modp.a = 1.2*1000/24;             % active transport rate (um/hr)
T05 = 51*24;                       % half-life (hours)
modp.g = 1/(T05/log(2));           % decay rate (per hour)
k1 = 0.2;                          % alpha_twid_h
modp.e0 = modp.g*modp.sg/(k1*modp.c0*modp.rg*modp.a);  % soma flux rate
% calculated parameters
[calcp] = CMNG_calcparams(simp, modp);
% run model for jmax time steps, linear ICs, no retraction
%[C1, C01, CN1, l1] = CMNG_rundyn(simp, modp, calcp, 0, 0);
%Ca1 = [C01 C1 CN1];
% get analytical steady-state values (note: not valid for ah=1)
[Cinfa1, linfa1] = CMNG_lCanal(simp, modp, calcp, 0);
linfa1 = linfa1*(modp.D/(modp.rg*modp.c0))/1000     % mm
Cinfa1 = Cinfa1*modp.c0;
ah1 = calcp.gamma*calcp.beta/(calcp.phi*calcp.alpha)
modp

% Run 3: Rabbit optic
modp.a = 2*1000/24;             % active transport rate (um/hr)
T05 = 14*24;                       % half-life (hours)
modp.g = 1/(T05/log(2));           % decay rate (per hour)
k1 = 0.55;                          % alpha_twid_h
modp.e0 = modp.g*modp.sg/(k1*modp.c0*modp.rg*modp.a);  % soma flux rate
% calculated parameters
[calcp] = CMNG_calcparams(simp, modp);
% run model for jmax time steps, linear ICs, no retraction
%[C1, C01, CN1, l1] = CMNG_rundyn(simp, modp, calcp, 0, 0);
%Ca1 = [C01 C1 CN1];
% get analytical steady-state values (note: not valid for ah=1)
[Cinfa1, linfa1] = CMNG_lCanal(simp, modp, calcp, 0);
linfa1 = linfa1*(modp.D/(modp.rg*modp.c0))/1000     % mm
Cinfa1 = Cinfa1*modp.c0;
ah1 = calcp.gamma*calcp.beta/(calcp.phi*calcp.alpha)
modp

% Run 4: Mouse optic
modp.a = 0.6*1000/24;             % active transport rate (um/hr)
T05 = 20*24;                       % half-life (hours)
modp.g = 1/(T05/log(2));           % decay rate (per hour)
k1 = 0.65;                          % alpha_twid_h
modp.e0 = modp.g*modp.sg/(k1*modp.c0*modp.rg*modp.a);  % soma flux rate
% calculated parameters
[calcp] = CMNG_calcparams(simp, modp);
% run model for jmax time steps, linear ICs, no retraction
%[C1, C01, CN1, l1] = CMNG_rundyn(simp, modp, calcp, 0, 0);
%Ca1 = [C01 C1 CN1];
% get analytical steady-state values (note: not valid for ah=1)
[Cinfa1, linfa1] = CMNG_lCanal(simp, modp, calcp, 0);
linfa1 = linfa1*(modp.D/(modp.rg*modp.c0))/1000     % mm
Cinfa1 = Cinfa1*modp.c0;
ah1 = calcp.gamma*calcp.beta/(calcp.phi*calcp.alpha)
modp

% Run 5: Goldfish optic
modp.a = 0.4*1000/24;             % active transport rate (um/hr)
T05 = 67*24;                       % half-life (hours)
modp.g = 1/(T05/log(2));           % decay rate (per hour)
k1 = 0.85;                          % alpha_twid_h
modp.e0 = modp.g*modp.sg/(k1*modp.c0*modp.rg*modp.a);  % soma flux rate
% calculated parameters
[calcp] = CMNG_calcparams(simp, modp);
% run model for jmax time steps, linear ICs, no retraction
%[C1, C01, CN1, l1] = CMNG_rundyn(simp, modp, calcp, 0, 0);
%Ca1 = [C01 C1 CN1];
% get analytical steady-state values (note: not valid for ah=1)
[Cinfa1, linfa1] = CMNG_lCanal(simp, modp, calcp, 0);
linfa1 = linfa1*(modp.D/(modp.rg*modp.c0))/1000     % mm
Cinfa1 = Cinfa1*modp.c0;
ah1 = calcp.gamma*calcp.beta/(calcp.phi*calcp.alpha)
modp

% Run 6: Blue whale (Alvarez et al)
modp.a = 7*1000/24;             % active transport rate (um/hr)
T05 = 730*24;                       % half-life (hours)
modp.g = 1/(T05/log(2));           % decay rate (per hour)
k1 = 0.14;                          % alpha_twid_h
modp.e0 = modp.g*modp.sg/(k1*modp.c0*modp.rg*modp.a);  % soma flux rate
% calculated parameters
[calcp] = CMNG_calcparams(simp, modp);
% run model for jmax time steps, linear ICs, no retraction
%[C1, C01, CN1, l1] = CMNG_rundyn(simp, modp, calcp, 0, 0);
%Ca1 = [C01 C1 CN1];
% get analytical steady-state values (note: not valid for ah=1)
[Cinfa1, linfa1] = CMNG_lCanal(simp, modp, calcp, 0);
linfa1 = linfa1*(modp.D/(modp.rg*modp.c0))/1000     % mm
Cinfa1 = Cinfa1*modp.c0;
ah1 = calcp.gamma*calcp.beta/(calcp.phi*calcp.alpha)
modp

% Run 7: Blue whale (Alvarez et al) - version 2: higher production
modp.a = 4*1000/24;             % active transport rate (um/hr)
T05 = 500*24;                       % half-life (hours)
modp.g = 1/(T05/log(2));           % decay rate (per hour)
k1 = 0.006;                          % alpha_twid_h
modp.e0 = modp.g*modp.sg/(k1*modp.c0*modp.rg*modp.a);  % soma flux rate
% calculated parameters
[calcp] = CMNG_calcparams(simp, modp);
% run model for jmax time steps, linear ICs, no retraction
%[C1, C01, CN1, l1] = CMNG_rundyn(simp, modp, calcp, 0, 0);
%Ca1 = [C01 C1 CN1];
% get analytical steady-state values (note: not valid for ah=1)
[Cinfa1, linfa1] = CMNG_lCanal(simp, modp, calcp, 0);
linfa1 = linfa1*(modp.D/(modp.rg*modp.c0))/1000     % mm
Cinfa1 = Cinfa1*modp.c0;
ah1 = calcp.gamma*calcp.beta/(calcp.phi*calcp.alpha)
modp

