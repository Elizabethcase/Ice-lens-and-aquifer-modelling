%given 

%function [U, Stefan, Pe, B] = converttoparams(

% constants from Tables 1 and 2
rhoI = 917; %ice density (ID in paper = rho); kg m-3
cp = 2050; %heat capacity of water; m2 s-2 K-1
L = 334000; %latent heat of water; m2 s-2
K = 2.1; %K-1
g = 9.806; %m s-2
gamma = .07; %N m-1
dp = 10e-4; %m
mu = 10e-3; %Pa s
Q0 = 200; %W m-2 or kg s-3
M = 6e-4; % kg m-2 s-1
l = 20.6; %m 
k0 = 5.6e-11; %m^2
t0 = 3.15e7; % number of seconds in a year
Sw = 292; %net shortwave radiation; W m-2
alpha = 0.6; %albedo
emis = 0.97; %emissivity
sbc = 5.7e-8; %Stefan-Boltzmann constant; W m-2 K-4
Lw = 279; %longwave radiation; W m-2
chi = 10.3; %turbulant transfer coeff; W m-2 K-1 
a0 = 9.5e-9; %accumulation rate m s-1
Ta = 273-4.1137; %average air temp; K %updated from juneauIceField_weather_v1.0/LVL2/juneauicefieldCamp18AWS_daily_LVL2.csv'
Tm = 273; %melting temp = 0 C

%shortwave from PoG 4th ed pg 143
th = 59.34211; %latitude of C18
del = ;%solar declination
h = ;%hour
cosz = sin(th)*sin(del)+cos(th)*cos(del)*cos(h);
Psy = Psi0^(P/P0/cosz);
Es = Es0 * cosz*Psy;

%relevant conversions
h = chi + 4*emis*sbc*Tm^3; %represents turbulant heat and outgoing longwave radiation; 14.8; %m-2 K-1 in paper
dT = Q0/h;

U = rhoI*g*k0*t0/l/mu; 
Stefan = L/cp/dT;
