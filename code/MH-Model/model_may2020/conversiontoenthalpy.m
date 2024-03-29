function [H,W] = conversiontoenthalpy(T,phi,S,Stefan)
% This function takes in the temperature, porosity, and saturation
% and puts out the enthalpy and total water/ice mass 

% [T,phi,S] = conversiontotemperature(H,W,Stefan)

global phi_min 

phi(phi<phi_min)=phi_min;
W = 1-phi+phi.*S;
H = W.*T+Stefan*phi.*S;

end
