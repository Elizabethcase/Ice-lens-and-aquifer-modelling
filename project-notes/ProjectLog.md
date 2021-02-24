	### Feb 23

	Created Simplified_RunCode_2020_10_20.m w/

```
no compaction
timestep is a set step, not a fraction of T
accumulation is normalized inside function
thetaOpt: optional 'seasonal' or 'constant' temperature cycle
accOpt: optional 'seasonal' or 'constant' acccumulation rate
phiOpt: initialized porosity profile ...
           'exponential', 
           'gausuni', (gauss fun + surface density)...
           'gausexp', (gauss fun + exponetial)
           'uniform', (=phi0), or ...
           'ice lens uni', (1 ice lens in uniform (=phi0) snowpack)...
           'ice lens exp', (1 ice lens in exponentially decreasing phi)
 R: choose vol water input from the top
 ```

Temp + acc: can choose from seasonal or uniform boundary condition
```
% Temperature
if strcmp('seasonal',thetaOpt)
    EbarFun = @(tau)Qbar-cos(2*pi*tau); %Q0 = 200
elseif strcmp('constant',thetaOpt)
    EbarFun = @(tau) Qbar;
end
% Accumulation
if strcmp('seasonal',accOpt)
    % Very simple seasonal accumulation (max in winter, min in summer)
    Abar = @(tau)AccumulationRate/2-AccumulationRate*sin(2*pi*tau); 
elseif strcmp('constant',accOpt)
    Abar = @(tau)AccumulationRate; % Accumulation
end
```

For porosity, I want to explore many options: initializing with an exponential, with a uniform profile, and then adding in an ice lens, either with a step change in porisity (ice lens..) or smoothly going to ice (gaussian function), per Colin's suggestion

```
% Porosity
if strcmp('exponential',phiOpt)
    Lphi = 5; %almost zero by z = 20 (~=ell, quite warm, equiv to juneau, also try increasing)
    phi = phi0 .* exp(-xgrid .* ell/Lphi);
elseif strcmp('gaussuni',phiOpt)
    aG = -0.54; %phi (close to ice density)
    bG = ell/2;
    cG = 1;
    phi = phi0 + aG * exp(-((xgrid .* ell - bG).^2)./(2*cG^2));
elseif strcmp('gaussexp',phiOpt)
    Lphi = 10; %almost zero by z = 20 (~=ell, quite warm, equiv to juneau, also try increasing)
    phi = phi0 .* exp(-xgrid .* ell/Lphi);
    %currently  just constant density + gaussian, but could also add
    %exponential instead
    aG = -phi(50)*.9; %phi (close to ice density)
    bG = ell/2;
    cG = 1;
    phi = phi + aG * exp(-((xgrid .* ell - bG).^2)./(2*cG^2));
elseif strcmp('uniform',phiOpt)
    phi = phi0 .* ones(N,1); 
elseif strcmp('ice lens uni',phiOpt)
    phi = phi0 .* ones(N,1);
    lensThickness = 2; % # spatial steps;
    zLoc = 5;
    [v,index] = min(abs(zLoc - xgrid*ell));
    phi(index:index + lensThickness) = 0;
elseif strcmp('ice lens exp',phiOpt)
    Lphi = 10; %almost zero by z = 20 (~=ell, quite warm, equiv to juneau, also try increasing)
    phi = phi0 .* exp(-xgrid .* ell/Lphi);
    lensThickness = 2; % # spatial steps;
    zLoc = 5;
    [v,index] = min(abs(zLoc - xgrid*ell));
    phi(index:index + lensThickness) = 0;
end 
```

Currently getting singular matrix error, probably due to how I'm initializing phi, will take another go tomorrow (Feb 24)