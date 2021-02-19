function [Fp,Fm] = SaturationDiffusiveFlux(Diffusivity,Porosity,Saturation,dt,dx,Constants)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% [Fp,Fm] = SaturationDiffusiveFlux(Diffusivity,Porosity,Saturation,dt,dx,Constants)
% Diffusive flux with input arguments: Diffusivity as a function handle of
% phi and S (two input arguments) and Porosity and
% Saturation, which come from the function
% [T,phi,S]=conversiontotemperature(H,W,Stefan).
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

Porositym = [Porosity(1); Porosity(1:(end-1))];
Porosityp = [Porosity(2:end); Porosity(end)];

Saturationm = [Saturation(1); Saturation(1:(end-1))];
Saturationp = [Saturation(2:end); Saturation(end)];

Dm = Diffusivity(Porositym,Saturationm);
D = Diffusivity(Porosity,Saturation);
Dp = Diffusivity(Porosityp,Saturationp);

fim = (2/dx)*(((1./Dm)+(1./D)).^(-1));
fip = (2/dx)*(((1./D)+(1./Dp)).^(-1));
Fm = Constants.*fim.*(Saturation-Saturationm);
Fp = Constants.*fip.*(Saturationp-Saturation);
end

