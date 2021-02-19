function [Fp,Fm] = EnthalpyDiffusiveFlux(Diffusivity,Temperature,Porosity,Saturation,dt,dx,Constants,Stefan)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% [Fp,Fm] = EnthalpyDiffusiveFlux(Diffusivity,Temperature,Porosity,Saturation,dt,dx,Constants,Stefan)
% Diffusive flux with input arguments: Diffusivity as a function handle of
% phi and S (two input arguments) and temperature, porosity, and
% Saturation, come from the function
% [T,phi,S]=conversiontotemperature(H,W,Stefan).
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Here the diffusivity includes the (T+Stefan) term.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% Temperaturem = [Temperature(end); Temperature(1:(end-1))];
% Temperaturep = [Temperature(2:end); Temperature(1)];
% Porositym = [Porosity(end); Porosity(1:(end-1))];
% Porosityp = [Porosity(2:end); Porosity(1)];
% Saturationm = [Saturation(end); Saturation(1:(end-1))];
% Saturationp = [Saturation(2:end); Saturation(1)];

Temperaturem = [Temperature(1); Temperature(1:(end-1))];
Temperaturep = [Temperature(2:end); Temperature(end)];

Porositym = [Porosity(1); Porosity(1:(end-1))];
Porosityp = [Porosity(2:end); Porosity(end)];

Saturationm = [Saturation(1); Saturation(1:(end-1))];
Saturationp = [Saturation(2:end); Saturation(end)];

Dm = Diffusivity(Porositym,Saturationm).*(Temperaturem+Stefan);
D = Diffusivity(Porosity,Saturation).*(Temperature+Stefan);
Dp = Diffusivity(Porosityp,Saturationp).*(Temperaturep+Stefan);

fim = (2/dx)*(((1./Dm)+(1./D)).^(-1));
fip = (2/dx)*(((1./D)+(1./Dp)).^(-1));
Fm = Constants.*fim.*(Saturation-Saturationm);
Fp = Constants.*fip.*(Saturationp-Saturation);
end

