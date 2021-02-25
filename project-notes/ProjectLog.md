### Feb 25

Drawing from MH paper, use experimental parameters & boundary conditions laid out in 3.1 (rainfall onto cold snow, p5)

- no accumulation
- no compaction
- S = 0 (dry snow at surface)
- phi = phi0
- T = T_inf < T_m

On the other hand, while they ensure porosity is large enough that the snow never saturates, I want the snow to saturate (and can force it to saturate at the induced gaussian ice lens which has a phi \~=0.1 (or less).

Just imposing a rain volume at the surface (e.g. constant influx of water), is not enough (no saturation anywhere, no evidence of surface melt input.) Need to look at ind. plots. Is the water showing up anywhere? If not, why not?

<img title='line plot' alt="Line plot of enthalpy, temp, sat, total water, porosity, water pressure, water flux" src=https://github.com/Elizabethcase/Ice-lens-and-aquifer-modelling/blob/main/project-notes/imgs/20210225/time_step_rain_input.png width=50%>

### Feb 24

I think the code is broken because of how I'm forcing the initialization of phi... I believe it needs to be done when setting W and H, which, at this point, just create a uniform phi of value phi0

```
W = (1-phi0)*ones(N,1); % take in W from above
H = W.*Qbar; % take in H from above
```

In Appendix B, 

$$ \mathcal{W} = 1 - \phi + \phi*S \text{ (total water)}$$ $$ \mathcal{H} = \mathcal{W}*\theta + \mathcal{S}*phi*S \text{ (enthalpy)}$$

where

$$ \mathcal{S} = \frac{\mathcal{L}}{c_p \Delta T} \text{ (Stefan number)}$$ $$ S = max\{0, \frac{\mathcal{H}}{\rho \mathcal{L} \phi}\} \text{ (saturation)} $$ $$ \phi = 1-\mathcal{W}+max\{0,\frac{\mathcal{H}}{\rho \mathcal{L}}\} \text{ (porosity)} $$

Tried changing phi0 to phi
      
```
W = (1-phi); % take in W from above
H = W.*Qbar; % take in H from above
```

But this did not solve the issue. Going back to double check that code still works w/ no phi initialization.

It does not.. what else have I broken? OK.

OK I had changed `dt` to be a timestep (instead of a fraction of the total time), but hadn't updated the rest of the code. That was breaking.

`W = (1-phi)` is still a more consistent way of initializing phi than forcing it with 

```
if n == 1
	phi_nm1 = phi;
end
```

And this seems to work! or at least run, and the normalized values look reasonable. Reaches steady state by ~ 10 yrs. (accumulation and temperature are constant, compaction is off, no water input, so porosity goes to phi0 as expected)

Next up, 
- [ ] thinking about boundary conditions at the surface with rain fall (e.g. a constant water input) 
- [ ] getting a feel for what happens as water is added at the surface with and without an ice lens
- [ ] initialize a temperature profile that's warm at the surface, cold deeper down.
- [ ] turn off accumulation?

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