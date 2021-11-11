function Shear = CompactionFunction(xgrid,W,phi,Theta,A,nglen,Accumulation,type)
% This function takes in W, phi, and theta and outputs the compaction rate
% xcelledges are the cell edge values
% phi is the snow porosity
% W is the total water, i.e. sum of ice and liquid water
% Theta is the snow temperature
% A is the snow viscosity softness
% nglen is the snow viscosity exponent
switch type
    case 'none'
        Shear = zeros(size(phi)); % dw/dz
        
    case 'poreclosure'
        IcePressure = cumtrapz(xgrid,W);
        Shear = 2.*A.*(phi./(1-phi)).*IcePressure.^nglen; % dw/dz
        
    case 'empirical'
        rhoi = 917; % ice density
        rhow = 1000; % water density
        f = 1-(550/rhoi); % transition porosity
        R = 8.314; % gas constant
        E0 = 10160; % activation energy
        E1 = 21400; % activation energy
        T = 273.15+(200/14.8)*Theta; % temperature in kelvin
        omega = 1/(3600*24*365); % 1/(seconds in a year)
        ell = 200/(omega*917*334000); % lengthscale
        k0 = 11; % H&L1980 constant
        k1 = 575; % H&L1980 constant
        A0 = (rhoi./rhow)*Accumulation*ell; % water equivalent accumulation per year % % ACCUMULATION TO MASS RATE
        A1 = (rhoi./rhow)*Accumulation*ell; % water equivalent accumulation per year % % ACCUMULATION TO MASS RATE
        a = 1; % H&L1980 constant
        b = 1/2; % H&L1980 constant
        c0 = (k0*(A0.^a)).*exp(-E0./(R*T));
        c1 = (k1*(A1.^b)).*exp(-E1./(R*T));
        C = c0.*(phi>=f)+c1.*(phi<f);
        Shear = C.*(phi./(1-phi)); % dw/dz
end
end

