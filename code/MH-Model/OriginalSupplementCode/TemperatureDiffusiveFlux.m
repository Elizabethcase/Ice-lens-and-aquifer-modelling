function [Fp,Fm] = TemperatureDiffusiveFlux(TotalWater,Temperature,dt,dx,Constants)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% [Fp,Fm] = TemperatureDiffusiveFlux(TotalWater,Temperature,dt,dx,Constants)
% Diffusive flux with input arguments: TotalWater which is W and 
% temperature, which comes from the function
% [T,phi,S]=conversiontotemperature(H,W,Stefan).
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

Temperaturem = [Temperature(1); Temperature(1:(end-1))];
Temperaturep = [Temperature(2:end); Temperature(end)];

TotalWaterm = [TotalWater(1); TotalWater(1:(end-1))];
TotalWaterp = [TotalWater(2:end); TotalWater(end)];

Dm = TotalWaterm;
D = TotalWater;
Dp = TotalWaterp;

fim = (2/dx)*(((1./Dm)+(1./D)).^(-1));
fip = (2/dx)*(((1./D)+(1./Dp)).^(-1));
Fm = -Constants.*fim.*(Temperature-Temperaturem);
Fp = -Constants.*fip.*(Temperaturep-Temperature);
end

