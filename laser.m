clc;
clear all;
close all;
 
%% input value
c=3e8;
lamda_0= 1200e-9;
v_0=c/lamda_0;
n2=3.3208; %active layer
n1=1;
n=n2;

%delv

delv=1.5e9 ;
%gth 
y=25;%loss coefficient( cm^-1)
L=2e-4;%length
w=8e-6;
d=0.2e-6;
A=w*d;
R1=((n2-n1)/(n2+n1))^2;
R2=R1;
gt=y+(1/(2*L*100))*log(1/(R1*R2));
gth=2*gt;
nth=1.5*10^24;

%tr
B=2*10^-16;
tr=(B*nth)^-1;

%Tsp
% Br=10^-11;
% n0=10^11;
% p0=10^8;
% del_n=3*10^18;
%tsp=(Br*(n0+p0))^-1;


delN=4e15;
%% gain coefficient
g_0=delN*((c*c)/(8*pi*v_0*v_0*n*n*tr*delv));

%tau_ph
Tph=n/(c*gth*100);

%nth
% nth=(8*pi*v_0*v_0*tr*delv*n*n*gth)/(c*c);

%Ith
e=1.6e-19;

Ith=(nth*e*L*w*d)/tr ;
Ith=Ith*1000;
 
%% guassian lineshape 
v_s=linspace((v_0-delv),(v_0+delv),1000);
v_m=mean(v_s);
v_std=std(v_s);
g_gauss= g_0*gaussmf(v_s,[v_std v_m]);
w_s=c./(v_s*1e-9);
plot((w_s),(g_gauss),'linewidth',2);
xlabel('wavelength (nm)')
ylabel ('gain coefficient')
hold on

v_min= v_0-delv/2;
v_max= v_0+delv/2;
w_max=c/v_min;
w_min=c/v_max;
m_max=round(2*n*L/w_min);
m_min=round(2*n*L/w_max);
m=m_min:1:m_max;
w=2*n*L./m;
w=w(find(w<=w_max));
w=w(find(w>=w_min));
w=w/1e-9;
g_discreet=interp1(w_s,g_gauss,w);
stem(w,g_discreet,'*','linewidth',3);

figure

%PowervsI

h=6.626e-34;
e=1.6e-19;
I = 0:1:100
P = zeros(length(I))
for i=1:length(I);

P(i)=((h*c*c*Tph*(1-R1))/(2*e*n*lamda_0*L))*(i-Ith);
end

plot(I,(P),'b','linewidth',2);
xlabel('current(mA)')
ylabel ('Power(mW)')

%I vs V

%% constants
%Na=1e15;
Na=4e16;
Nd=5e16;
%Nd=5e17;
ni=4.4*10^9;
npo=ni^2/Na;
pno=ni^2/Nd;
deln=1.5e18;
un=3900;
up=90;
K=1.38e-23;
T=300;
q=1.6e-19;
Dn=(K*T*un)/q;
Dp=(K*T*up)/q;
Br=1e-11;
tn=1/(Br*(Nd+pno+deln));
tp=1/(Br*(npo+Na+deln));
Ln=(Dn*tn)^0.5*1e-2;
Lp=(Dp*tp)^0.5*1e-2;
nf=1;
V=0:.001:1.4;
for i=1:length(V)
    I(i)=q*A*((Dn*npo/Ln)+(Dp*pno/Lp))*exp(((q*V(i))./(nf*K*T))-1);
end
plot(V,(I),'linewidth',2);
xlabel('Voltage (V)');
ylabel('Current (mA)');

title('I-V Characteristics of Laser');