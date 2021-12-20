function [T,phi,S] = conversiontotemperature(H,W,Stefan)
% This function takes in the enthalpy and total water/ice mass and spits
% out the temperature, porosity, and saturation
% [T,phi,S] = conversiontotemperature(H,W,Stefan)

T = min(0,(H./W)); 
    % max temp = 0 

%phi = max(1-W+(1/Stefan)*max(0,H),0);
phi = 1-W+(1/Stefan)*max(0,H);
phi(phi<1e-4)=1e-4;

% Sc = max(0,(H./(Stefan*phi)));
% S = max(0,H)./(max(0,H)+Stefan*(1-W));
S = max(min(1,max(0,H)./(max(0,H)+Stefan*(1-W))),0);
%solve for S for phi = phi_min

%if phi == phi_min
%    S(H,W) = 

    % inner min -- 
        % 1-W = phi - phi*S --> amount of unfilled pore space
        % 1-W should never be > 1; 0 if pore space saturated
        % Stefan = 12, so when (1-W) > 1/12, min == 1
    % outer max --
        % pore space is either saturated or it isn't
        % if min() = 1, saturated 
        % if min() = 0, unsaturated

end