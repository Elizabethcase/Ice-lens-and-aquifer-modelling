% non-normalized Q(t) given T

Ta = -4 + 273;
%Q = QgivenT(Ta)
Q0 = 200;

Qnorm = Q/Q0

%function Q=QgivenT(Ta)
    alpha = 0.6; % albedo
    Sw = 292; % net shortwave radiation (W m-2)
    epsilon = 0.97; %emissivity
    sig = 5.7 * 10^-8; %stefan-boltzman (W m-2 K-4)
    Lw = 279; %Longwave radiation (W m-2)
    chi = 10.3; % Turbulant transfer (W m-2 K-1)
    rhow = 1000;
    R = 8.05 * 10^-10; % rainfall rate? 8.05 * 10^-10 = 1 inch/year (true value uniknown
    a = 9.26 * 10^-8; % accumulation rate (m.w.e s-1) 9.26 * 10^-8;  = 3 m.w.e./yr
    ci = 2108; % heat capacity of ice (J kg-1 K-1)
    cw = 4179; % heat capacity of water (J kg-1 K-1)
    Tm = 273; % melting temperature
    h = 14.8; %W m-2 K-1 (chi + 4*epsilon*sig*Tm^4
    
    Q = (1-alpha) * Sw + Lw - epsilon * sig * Tm^4 ...
        + chi * (Ta - Tm) + rhow * ci * a * (Ta-Tm) ...
        + rhow * cw * R * (Ta-Tm);

    tau = 
    Qbar = Q - Q0*cos(2*pi*tau);
%end