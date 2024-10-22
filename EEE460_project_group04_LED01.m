clc
clear all
format long

set(0,'DefaultAxesFontName', 'Times');
set(0,'DefaultAxesFontSize', 11);


%% Constant 
m0     = 9.1e-31;
q      = 1.6e-19;
kb     = 1.38e-23;
h      = 6.6626e-34;
h_cut  = h/(2*pi);
c      = 3e8;


%% Simulation Parameter Input
lambda = 410:5:700;          % nm
T      = 300;                 % kelvin

kbT    = kb*T/q;              % eV
eV     = 1243./lambda;        % eV


%% Material Properties((AlxGa1-x)0.5In0.5P) x = 0.36
me          = 0.094*m0;               % elec effective mass
mh          = 0.62*m0;                % hole effective mass
mu_e        = 100;                    % electron mobility
mu_h        = 7;                      % hole mobility
Eg          = 2.13;                   % Bandgap
ni          = 2*(2*pi*kb*T/h^2)^1.5*(me*mh)^0.75*exp(-Eg/(2*kbT));                   
%ni         = 4.5e17                  % Intrinsic carrier conc.
Br          = 1e-10;                  % Coefficient for band-band recombination
sr          = 1e-14;                  % conc. of recombination centers 2
Nt          = 1e2;                    % capture cross-section
nr1         = 2.5;                    % refractive index of semiconductor
nr2         = 1;                      % refractive index of air
nf          = 1;                      % ideality factor
A           = 1e-2;                   % Area in cm^2
Nd          = 5e17;                   % n type doping profile
Na          = 1e15;                   % p type doping profile
del_n       = 1e18;                   % excess carrier in p side

del_p       = del_n;                  % excess carrier in n side
p_No        = ni^2/Nd;                % hole conc. in n type
n_Po        = ni^2/Na;                % elec conc. in p type
Efn         = -kbT*log(Nd/ni);        % fermi-level for elec
Efp         = kbT*log(Na/ni);         % fermi-level for hole
Dn          = kbT*mu_e;               % diffusion coeff of elec 
Dp          = kbT*mu_h;               % diffusion coeff of hole
tau_r_N     = 1/(Br*(Nd+p_No+del_p)); % diffusion lifetime of elec
tau_r_P     = 1/(Br*(n_Po+Na+del_n)); % diffusion lifetime of hole
Ln          = (Dn*tau_r_N)^0.5*1e-2;  % diffusion length of elec
Lp          = (Dp*tau_r_P)^0.5*1e-2;  % diffusion length of hole
mr          = (me*mh)/(me+mh);        % reduced mass
vth         = sqrt((3*kbT*q)/mr)*1e2; % thermal velocity of elec
tnr         = (sr*vth*Nt)^-1;         % indirect recombination lifetime


%% Spontaneous Emmission Rate or Injection Electroluminiscent Rate per unit volume
a1     = (((2*mr)^(3/2))/(2*pi^2*h_cut^3*tau_r_P)).*(eV-Eg).^(0.5);
a2     = exp((Efn-Efp-Eg)/kbT);
a3     = exp(-(eV-Eg)./kbT);
rsp    = a1.*a2.*a3;

figure;
plot(lambda, rsp*1e-6*1e18, 'linewidth', 2)
xlabel('\lambda(nm)')
ylabel('r_{sp} (\times 10^{-18} (s^{-1}.(eV)^{-1}.cm^{-3}))');
title('spontaneous emission rate vs wavelength ((Al_{0.36}Ga_{0.64})_{0.5}In_{0.5}P)', 'FontSize', 10);
grid on;
grid(gca,'minor');

 
%% Peak Wavelength & Linewidth
rsp1      = rsp;
rsp       = real(rsp);
halfmax   = (max(rsp)+min(rsp))/2;
index1    = find(rsp >= halfmax, 1, 'first');
index2    = find(rsp >= halfmax, 1, 'last');
linewidth = lambda(index2) - lambda(index1);
lambda0   = lambda(find(rsp == max(rsp)))*1e-9;


%% Injection Efficiency
nin         = (Dn*n_Po/Ln)/((Dn*n_Po/Ln)+(Dp*p_No/Lp)) 

%% Radiative Recombination Efficiency(IQE)
nr          = 1/(1+(tau_r_N/tnr)) 

%%  Extraction Efficiency(EE)
crit_ang    = asin(nr2/nr1);
Ft          = (1/4)*(nr2/nr1)^2*(1-((nr1-nr2)/(nr1+nr2))^2);
ne          = (1/2)*(1-cos(crit_ang))*Ft 

%% File Input for Luminous Efficiency
C           = dlmread('sensitivity_GaAs.txt'); %sensitivity file
lambda      = C(:,1);
sensitivity = C(:,2);
emission    = rsp;

%% Luminous Efficiency
v1          = emission*sensitivity;
V           = sum(v1);
P           = sum(emission);
nl          = V/P


%% Output Optical Power(L) as a Function of Forward Current(I)               
hv          = h*c/lambda0;
n0          = nin*nr*ne*nl;

I           = 1:5:50;
P0          = (n0*hv/q)*I;

figure
plot(I, P0*1e3, 'linewidth', 2)
grid on;
xlabel('forward current, I(mA)')
ylabel('output optical power, L(mW)')
title('Light-Current Characteristics Of LED((Al_{0.36}Ga_{0.64})_{0.5}In_{0.5}P)');
grid(gca,'minor')


%% Determination of Current
V           = 0:0.01:2;

% reverse saturation current
Is          = q*A*((Dn*n_Po)/Ln + (Dp*p_No)/Lp);
I           = Is.*exp((V)./(nf*kbT))*1e3;    % current in mA

figure;
plot(V, I,'LineWidth', 2);
xlabel('Applied Forward Bias Voltage, V');
ylabel('I(mA)');
title('I-V Characteristics Of a n^{+}-p LED ((Al_{0.36}Ga_{0.64})_{0.5}In_{0.5}P)');
%xlim([0 3]);