function [Zu, Zl] = solutionto33(phi, phi0, R, Z0, Z1,xgrid)

beta = 2;
rhoI = 917;
mu = 10e-3;
k0 = 5.6e-11;
g = 9.81;

Zu_prev = Z1;
Zl_prev = Z1;


[val, Iu] = min(abs(xgrid-Zu));
[val, Il] = min(abs(xgrid-Zl)); 

qs = 3*k0*(phi0^3)*rhoI*g*(Zu_prev - Zl_prev)./...
        mu*Z0*(exp(3*Zu/Z0) - exp(3*Zl/Z0));

Zu = Zu_prev + dt.*(...
        (qs-R)./...
        (phiU.*(1-((mu * R * exp(3*Zu/Z0))/(rhoI*g*k0*phi0^3))^(1/beta))));

Zl = qs./phil;

Zu = Zu_prev;
Zl = Zl_prev;

end