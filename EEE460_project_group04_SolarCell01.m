clc
clear all
format long

set(0,'DefaultAxesFontName', 'Times');
set(0,'DefaultAxesFontSize', 11);


%% Constant 
m0      = 9.1e-31;
q       = 1.6e-19;
c       = 3e8;
kb      = 1.38e-23;
h       = 6.6626e-34;
h_cut   = h/(2*pi);
T       = 300;                                  % Temperature(K)
kbT     = kb*T;                                 % (eV)

%% Material Properties((Si) 
me          = 0.98*m0;                % elec effective mass
mh          = 0.49*m0;                % hole effective mass
mu_e        = 1400;                   % electron mobility
mu_h        = 450;                    % hole mobility
Eg          = 1.12;                   % Bandgap
%ni         = 2*(2*pi*kb*T/h^2)^1.5*(me*mh)^0.75*exp(-Eg/(2*kbT));                   
ni          = 1.5e10;                 % Intrinsic carrier conc.
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


%% Solar Cell Properties(Si)
Rs      = [0.01 0.02 0.03 0.04 0.05];           % Series Resistance(Ohm)
V       = 0: 0.001: 1.2;                        % Varing Voltage(V)
n       = 1;                                    % Ideality Factor: Ideal value = 1
I       = 1000;                                  % Radiation Intensity(W/m^2)
Iph     = 20;                                   % PhotoCurrent(mA)
Is      = q*A*((Dn*n_Po)/Ln + (Dp*p_No)/Lp);    % Reverse Saturation Current(mA)


%% Calculation
for i=1:length(Rs)
    for j=1:length(V)
        % Diode Current(mA as Iph & Io both in mA)
        f            = @(i_t) - i_t - Iph + Is*(exp((q*(V(j) - i_t*Rs(i))/(n*kbT))) - 1); 
        i_t          = fzero(f, 20);           % solving non-linear equation wth initial value Iph
        I_total(i,j) = i_t;
    end
end

P       = (-I_total).*V;                        % Power(mW = mA*V)

[x, i]  = min(abs(I_total), [], 2);             % Find indices where I_total = 0
V_oc    = V(i)';                                % Open Circuit Voltage, V_oc = V(I_tolal = OA)
I_sc    = 10*ones(length(Rs), 1);               % Short Circuit Current, I_sc = I_total(0V) = -Iph + 0
MaxP    = (max(P'))';
FF      = MaxP./(V_oc.*I_sc);                   % Fill Factor

figure
plot(V, I_total, 'linewidth', 2);
hold on;
grid on;
grid(gca,'minor');
xlabel('Voltage, V(V)');
ylabel('Current, I_{total}(mA)');
title('IV Characteristics(Si Solar Cell)');
legend('R_s = 0.01\Omega','R_s = 0.02\Omega','R_s = 0.03\Omega','R_s = 0.04\Omega','R_s = 0.05\Omega', 'Location', 'best');
xlim([0 max(V_oc)]);
ylim([-20 1]);
yline(0);

figure
plot(V, P, 'linewidth', 2);
hold on;
grid on;
grid(gca,'minor');
xlabel('Voltage, V(V)');
ylabel('Power, P(mW)');
title('PV Characteristics(Si Solar Cell)');

legend('R_s = 0.01\Omega','R_s = 0.02\Omega','R_s = 0.03\Omega','R_s = 0.04\Omega','R_s = 0.05\Omega', 'Location', 'best');
xlim([0 max(V_oc)]);
ylim([0 max(MaxP)]);
yline(0);

figure
plot(Rs,  MaxP, 'linewidth', 2);
hold on;
grid on;
grid(gca,'minor');
xlabel('Series resistance, R_s(\Omega)');
ylabel('Maximum Power, P(mW)');
title('P_{max} vs R_s Graph(Si Solar Cell)');

figure
plot(Rs, V_oc, 'linewidth', 2);
hold on;
grid on;
grid(gca,'minor');
xlabel('Series resistance, R_s(\Omega)');
ylabel('Open Circuit Voltage, V_{oc}(V)');
title('V_{oc} vs R_s Graph(S Solar Cell)');

figure
plot(Rs, I_sc, 'linewidth', 2);
hold on;
grid on;
grid(gca,'minor');
xlabel('Series resistance, R_s(\Omega)');
ylabel('Short Circuit Current, I_{sc}(mA)');
title('I_{sc} vs R_s Graph(Ideal Solar Cell)');

figure
plot(Rs, FF, 'linewidth', 2);
hold on;
grid on;
grid(gca,'minor');
xlabel('Series resistance, R_s(\Omega)');
ylabel('Fill Factor, FF');
title('FF vs R_s Graph(Ideal Solar Cell)');