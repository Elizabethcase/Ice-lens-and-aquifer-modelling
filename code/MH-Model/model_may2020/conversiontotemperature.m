function [T,phi,S] = conversiontotemperature(H,W,Stefan)
% This function takes in the enthalpy and total water/ice mass and spits
% out the temperature, porosity, and saturation
% [T,phi,S] = conversiontotemperature(H,W,Stefan)
global phi_min
global Stefan

T = min(0,(H./W)); 
    % max temp = 0 

%phi = max(1-W+(1/Stefan)*max(0,H),0);
phi = 1-W+(1/Stefan)*max(0,H);
phi(phi<phi_min)=phi_min;

W_sat = W;
% W_sat(W_sat>1)=1;
H_sat = H;
% H_sat(phi==0) = T(phi==0);

% Sc = max(0,(H./(Stefan*phi)));
% S = max(0,H)./(max(0,H)+Stefan*(1-W));
S = max(min(1,max(0,H_sat)./(max(0,H_sat)+Stefan*(1-W_sat))),0);
%S(phi==phi_min) = H(phi==phi_min)./(Stefan*phi_min);

%W2 = 1-phi+phi.*S;
%H2 = W.*T+Stefan*phi.*S;
 
%solve for S for phi = phi_min

    % inner min -- 
        % 1-W = phi - phi*S --> amount of unfilled pore space
        % 1-W should never be > 1; 0 if pore space saturated
        % Stefan = 12, so when (1-W) > 1/12, min == 1
    % outer max --
        % pore space is either saturated or it isn't
        % if min() = 1, saturated 
        % if min() = 0, unsaturated

end