clc;
clear all;
close all;

%%
X = dlmread('InGaAs.txt');
lambda = X(:,1);  %um
n = X(:,2);
k = X(:,3);

alpha = 4*pi.*k./lambda;  %absorption coefficient 
r1 = (n-1).^2 + k.^2;
r2 = (n+1).^2 + k.^2;
R  = r1./r2; % reflectivity 

wp = 0.005; %um
wi = 1;  %um
ne = (1-R).*(exp(-1*alpha.*wp) - exp(-1*alpha.*(wp+wi)));  % quantum efficiency 
A  = 20;       
e  = 1.6e-19;
I0 = 500;   
h  = 6.626e-34;
c  = 3e8;
lambda = lambda.*1
Iph= ne.*A.*e.*I0.*lambda./(h*c);  %7.30(a) equation by Willson
figure(1)
plot(lambda,Iph)
xlabel('Wavelength(um)','FontSize',16);
ylabel('Photocurrent (Iph)','FontSize',16);
xline(1.2);
grid on
title(' Photocurrent (Iph) vs. wavelength ');
figure(2)
plot(lambda,ne.*lambda,'linewidth',2);
xlabel('Wavelength(um)','FontSize',16);
ylabel('Responsivity (R)','FontSize',16);
xline(1.2);
grid on
title(' Responsivity (R) vs. wavelength ');

RL = 10;
Vext = Iph.*RL;
figure(3)
plot(Iph,Vext,'linewidth',2);
xlabel('V (external)','FontSize',16);
ylabel('Current','FontSize',16);
grid on
title(' Voltage vs Current ');
%%Gain calculation
rate_of_e_flow = Iph./e;
w = 10; %um
l = 1; %um
rate_of_e_gen = w.*l.*ne.*I0.*lambda./(h*c);
G = rate_of_e_flow./rate_of_e_gen;
figure(4)
plot (lambda, G);

%% Constant 
m0     = 9.1e-31;
q      = 1.6e-19;
kb     = 1.38e-23;
h      = 6.6626e-34;
h_cut  = h/(2*pi);


%% Simulation Parameter Input
%lambda = 300:1:1000;          % nm
T      = 300;                 % kelvin
kbT    = kb*T/q;              % eV
eV     = 1243./lambda;        % eV


%% Material Properties(GaAsP)
me          = 0.05216*m0;               % elec effective mass
mh          = 0.164*m0;               % hole effective mass
mu_e        = 14775;                   % electron mobility  in cm^2v^-1s^-1
mu_h        = 416.26;                    % hole mobility  in cm^2v^-1s^-1
Eg          = 1.13;                   % Bandgap
ni          = 2.98e14;                 % Intrinsic carrier conc.
Br          = 1e-10;                  % Coefficient for band-band recombination
sr          = 1e-14;                  % conc. of recombination centers 
Nt          = 1e2;                    % capture cross-section
%nr1         = ;                    % refractive index of semiconductor
%nr2         = 1;                      % refractive index of air
nf          = 2;                    % ideality factor
A           = 20;                   % Area in cm^2
Nd          = 5e17;                   % n type doping profile
Na          = 1e15;                   % p type doping profile
del_n       = 1e18;                   % excess carrier in p side

del_p       = del_n;                  % excess carrier in n side
p_No        = ni^2/Nd;                % hole conc. in n type
n_Po        = ni^2/Na;                % elec conc. in p type
Dn          = kbT*mu_e;               % diffusion coeff of elec 
Dp          = kbT*mu_h;               % diffusion coeff of hole
tau_r_N     = 1/(Br*(Nd+p_No+del_p)); % diffusion lifetime of elec
tau_r_P     = 1/(Br*(n_Po+Na+del_n)); % diffusion lifetime of hole
Ln          = (Dn*tau_r_N)^0.5*1e-2;  % diffusion length of elec
Lp          = (Dp*tau_r_P)^0.5*1e-2;  % diffusion length of hole
mr          = (me*mh)/(me+mh);        % reduced mass
vth         = sqrt((3*kbT*q)/mr)*1e2; % thermal velocity of elec
tnr         = (sr*vth*Nt)^-1;         % indirect recombination lifetime

% reverse saturation current

Is          = q*A*((Dn*n_Po)/Ln + (Dp*p_No)/Lp);
%% Determination of Current
Vtot=-10:.1:0;
Iph1=4.6572e-3; %at 1100nm
Itot=-Iph1+Is*(exp(Vtot/kbT)-1);
figure(6)
plot(Vtot,Itot,'linewidth',2);
xlabel('Voltage','FontSize',16);
ylabel('Current','FontSize',16);
grid on
title(' Voltage vs Current ');

I_dark = Is;
% NEP =

Responsivity  = ne.*lambda;
NEP = (1./Responsivity) .* sqrt(2*e*(I_dark + Iph));
in = 2*e*sqrt((I_dark + Iph));
Detectivity = sqrt(A)./ NEP;