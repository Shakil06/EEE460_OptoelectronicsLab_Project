clc
clear all
format long


%% Constant 
m0          = 9.1e-31;
q           = 1.6e-19;
c           = 3e8;
kb          = 1.38e-23;
h           = 6.6626e-34;
h_cut       = h/(2*pi);

%% Simulation Parameter Input
%lambda     = 410:5:700;              % nm
T           = 300;                    % kelvin

kbT         = kb*T/q;                 % eV
%eV         = 1243./lambda;           % eV

%% Material Properties(In1-xGaxAs) x = 0.79
me          = 0.054*m0;               % elec effective mass
mh          = 0.4*m0;                 % hole effective mass
mu_e        = 6953;                   % electron mobility
mu_h        = 400;                    % hole mobility
Eg          = 1.13;                   % Bandgap
ni          = 6.3e11;                 % Intrinsic carrier conc.
Br          = 2e-16;                  % Coefficient for band-band recombination
sr          = 1e-14;                  % conc. of recombination centers 2
Nt          = 1e2;                    % capture cross-section
nr1         = 3.38;                   % refractive index of semiconductor
nr2         = 1;                      % refractive index of air
nf          = 1;                      % ideality factor
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

dN          = 4e15;
lambda0     = 1100e-9;                 
v_0         = c/lambda0;
delv        = 1.5e9;
L           = 0.5;
w           = 8e-6;
d           = 0.2e-6;
A           = w*d;
y           = 0.17;
Ci          = 2e-16;

%% R1=100%, R2=90%
R1          = 1;
R2          = 0.7;
gth         = y+(1/(2*L))*log(1/(R1*R2));

%% gain coefficient
Tph        = nr1/(c*gth*100);
nth        =  1/(Ci*Tph);
tr         = (Br*nth)^-1;
Ith        = 1000*(nth*q*L*A)/tr;

g_0        = dN.*((nr1.*c.*c)./(8*pi*v_0*v_0.*tr.*delv));

%% guassian lineshape 
v_s        = linspace((v_0-delv),(v_0+delv),1000);
g_gauss    = g_0'.*gaussmf(v_s,[std(v_s) v_0]);
w_s        = c./(v_s*1e-9);

figure;
plot(w_s, g_gauss', 'linewidth', 2);
title('g(\nu_0) vs \nu_0 graph(In_{0.21}Ga_{0.79}As)');
hold on;

%% cavity modes
w_max      = c/(v_0-delv*1/2);
w_min      = c/(v_0+delv*1/2);
m_max      = round(2.*L./w_min);
m_min      = round(2.*L./w_max);
m          = m_min:1:m_max;
w          = 2.*L'./m;
w          = w(find(w<=w_max & w>=w_min));
w          = w./1e-9;
g_discreet = interp1(w_s, g_gauss', w);

stem(w, g_discreet, '*', 'linewidth', 3);
xlabel('wavelength, \lambda(nm)');
ylabel ('gain coefficient, g(\nu_0)');
legend(sprintf('mode: %.0f', length(w)));
hold on;
yline(gth);

%% Determination of Current
V           = 0:0.001:1.3;

% reverse saturation current
Is          = q*A*((Dn*n_Po)/Ln + (Dp*p_No)/Lp);
I           = Is.*exp((V)./(nf*kbT))*1e3;    % current in mA

figure;
plot(V, I,'LineWidth', 2);
xlabel('Applied Forward Bias Voltage, V');
ylabel('I(mA)');
title('I-V Characteristics Of a QW LASER (In_{0.21}Ga_{0.79}As)');

%% Determination of Power
P =((h*c*c*Tph*(1-R2))/(2*q*nr1*lambda0*L)).*(I-Ith);
P(P<0) = 0;

figure;
plot(I, P, 'b', 'linewidth', 2);
xlabel('current(mA)')
ylabel ('Power(mW)')
title('P-I Characteristics Of a QW LASER (In_{0.21}Ga_{0.79}As)');
%% Slope Efficiency
n_slope     = P(1300)/(I(1300)-Ith); 

%% External Quantum Efficiency(EQE)
n_eqe       = q*P(1300)/(I(1300)*Eg); 

%% External Differential Quantum Efficiency(EDQE)
n_edqe       = q*P(1300)/((I(1300)-Ith)*Eg); 

%% Internal Quantum Efficiency(IQE)
n_iqe       = 1/(1+(tau_r_N/tnr)); 

%% Extraction Efficiency(EE)
n_ee          = (gth - y)/gth;

%% External Quantum Efficiency(PCE)
n_pce       = P(1300)/(I(1300)*V(1300)); 

%title('PI Characteristics(InGaAs LASER)', sprintf('eta_{slope} = %.4f, eta_{EQE} = %.4f, eta_{EDQE} = %.4f, eta_{IQE} = %.4f, eta_{EE} = %.4f, eta_{PCE} = %.4f', n_slope, n_eqe, n_edqe, n_iqe, n_ee, n_pce));