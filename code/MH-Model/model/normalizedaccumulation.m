% Takes in meters of snow, ice, or water per year, and outputs normalized
% rate based on given parameters
% If norm is given as units, returns in meters of snow per year

% a = meters of accumulation in given units
% units = m, mie, mwe, norm
% Q0 = 200; %energy forcing %W m-2 or kg s-3
% t0 = 3.15*10^7 % seconds per year s
% rhoI = 917 % ice denisty kg m-3
% L = 334000 %latent heat of water; m2 s-2
% phi0 = 0.64; % surface percent of void space

function AN = normalizedaccumulation(a, units, Q0, t0, rhoI, L, phi0)

l = Q0*t0/(rhoI*L);

if strcmp(units,'m')
    
    AN = a./l;
    
elseif strcmp(units, 'mie')
    
    AN = a./(l*(1-phi0));
    
elseif strcmp(units, 'mwe')
    
    AN = (a./(l*(1-phi0)))*1000/917;
    
elseif strcmp(units,'norm')
    
    AN = a.*l;
else
    
    disp(['Unit ' units ' unavailable as an input.']);
    
end

end