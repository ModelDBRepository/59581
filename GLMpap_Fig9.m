% Continuum Model for Neurite Outgrowth with Autoregulation
% Graham, Lauchlan & McLean Figure 9
% Degenerate elongation (no decay) with and without tubulin autoregulation
% Version 1.0 (BPG & DRM 7-2-05)

% Parameters

% simulation
simp.dt = 0.01;                  % time step
simp.tmax = 5000;                % simulation time
simp.datat = 1000;                % data collection time step
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
k1 = 0.5;
k2 = 0.00001;
modp.e0 = modp.g*modp.sg/(k1*modp.c0*modp.rg*modp.a);  % soma flux-source rate
theta = 0;                      % fractional autoregulation
modp.er = theta*modp.e0;        % soma tubulin autoregulation
modp.rdt = 0;                   % autoregulation time delay
modp.el = k2*modp.rg;           % growth cone flux-sink rate
modp.zl = k2*modp.sg;           % growth cone flux-source rate
modp.g = 0;                     % remove decay

% plot parameters
tfs = 12;   % title font size


doruns = 1;     % flag to run simulations
if (doruns == 1)
    
% Run 1: standard D & a; no autoregulation
% calculated parameters
[calcp] = CMNG_calcparams(simp, modp);
% run model for jmax time steps, linear ICs
[Cgs, C0gs, CNgs, lgs] = CMNG_run(simp, modp, calcp, -1, modp);
[tgs, Cgs, C0gs, CNgs, lgs] = CMNG_dimen(simp, modp, Cgs, C0gs, CNgs, lgs);  % dimensionalise
Cags = [C0gs Cgs CNgs];
    
% Run 2: set a to 50%
modp.a = 50;
% calculated parameters
[calcp] = CMNG_calcparams(simp, modp);
% run model for jmax time steps, linear ICs
[Cga, C0ga, CNga, lga] = CMNG_run(simp, modp, calcp, -1, modp);
[tga, Cga, C0ga, CNga, lga] = CMNG_dimen(simp, modp, Cga, C0ga, CNga, lga);  % dimensionalise
Caga = [C0ga Cga CNga];
    
% Run 3: set D to 50%
modp.a = 100;
modp.D = 15000;
% calculated parameters
[calcp] = CMNG_calcparams(simp, modp);
% run model for jmax time steps, linear ICs
[CgD, C0gD, CNgD, lgD] = CMNG_run(simp, modp, calcp, -1, modp);
[tgD, CgD, C0gD, CNgD, lgD] = CMNG_dimen(simp, modp, CgD, C0gD, CNgD, lgD);  % dimensionalise
CagD = [C0gD CgD CNgD];

% Run 4: standard D & a; with autoregulation
modp.D = 30000;
k1 = 0.4;
modp.e0 = 0.002*modp.sg/(k1*modp.c0*modp.rg*modp.a);  % soma flux-source rate
theta = 0.2;                    % fractional autoregulation
modp.er = theta*modp.e0;        % soma tubulin autoregulation
% calculated parameters
[calcp] = CMNG_calcparams(simp, modp);
% run model for jmax time steps, linear ICs, change to large growth
[Cgr, C0gr, CNgr, lgr] = CMNG_run(simp, modp, calcp, -1, modp);
[tgr, Cgr, C0gr, CNgr, lgr] = CMNG_dimen(simp, modp, Cgr, C0gr, CNgr, lgr);  % dimensionalise
Cagr = [C0gr Cgr CNgr];

end

% Plot results

subplot(3,3,1);
plot(tgs,lgs,'k-');
hold on;
plot(tga,lga,'k-.');
title('a','FontSize',tfs);
ylabel('Length (\mum)');
%legend('a=100','a=50',4);
subplot(3,3,4);
plot(tgs,C0gs,'k-');
hold on;
plot(tga,C0ga,'k-.');
ylabel('Soma Concentration (\muM)');
subplot(3,3,7);
plot(tgs,CNgs,'k-');
hold on;
plot(tga,CNga,'k-.');
xlabel('Time (hours)');
ylabel('Terminal Concentration (\muM)');
legend('a=100','a=50',4);

subplot(3,3,2);
plot(tga,lgs,'k-');
hold on;
plot(tgD,lgD,'k-.');
title('D','FontSize',tfs);
ylabel('Length (\mum)');
%legend('D=30k','D=15k',4);
subplot(3,3,5);
plot(tga,C0gs,'k-');
hold on;
plot(tgD,C0gD,'k-.');
subplot(3,3,8);
plot(tgs,CNgs,'k-');
hold on;
plot(tgD,CNgD,'k-.');
xlabel('Time (hours)');
legend('D=30000','D=15000',4);

subplot(3,3,3);
plot(tgs,lgs,'k-');
hold on;
plot(tgr,lgr,'k-.');
title('AR','FontSize',tfs);
%legend('no AR','AR',4);
subplot(3,3,6);
plot(tgs,C0gs,'k-');
hold on;
plot(tgr,C0gr,'k-.');
subplot(3,3,9);
plot(tgs,CNgs,'k-');
hold on;
plot(tgr,CNgr,'k-.');
xlabel('Time (hours)');
legend('no AR','AR',4);



