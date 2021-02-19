function [T,phi,S] = conversiontotemperature(H,W,Stefan)
% This function takes in the enthalpy and total water/ice mass and spits
% out the temperature, porosity, and saturation
% [T,phi,S] = conversiontotemperature(H,W,Stefan)
T = min(0,(H./W));
phi = max(1-W+(1/Stefan)*max(0,H),0);
% Sc = max(0,(H./(Stefan*phi)));
% S = max(0,H)./(max(0,H)+Stefan*(1-W));
S = max(min(1,max(0,H)./(max(0,H)+Stefan*(1-W))),0);
end