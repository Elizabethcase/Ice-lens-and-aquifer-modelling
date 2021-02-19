function [qp,qm,pressure] = FullySaturatedWaterPressure(U,k,pressure,dx,xgrid,N,A,phi,W,Theta,nglen,I,Accumulation,type)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Compute water pressure when the snowpack is fully saturated. I is an
% indicator function that determines whether to solve for the water pressure
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

phimin = 10^(-4);
I(phi<phimin)=0;

phim = [phi(1); phi(1:(end-1))]; % cell values to the left of cell i
phip = [phi(2:end); phi(end)]; % cell values to the right of cell i

Dm = U*k(phim);
D = U*k(phi);
Dp = U*k(phip);

fim = (2/dx)*(((1./Dm)+(1./D)).^(-1));
fip = (2/dx)*(((1./D)+(1./Dp)).^(-1));
fip(N)=0; % no pressure gradient
% full difference
S = [fim(2:N)', -(fim+fip)',fip(1:(N-1))'];
Ind = [2:N, 1:N, 1:(N-1)];
Jnd = [1:(N-1), 1:N, 2:N];
M = sparse(Ind,Jnd,S,N,N);
% positive flux
S = [zeros(1,N-1),-fip',fip(1:(N-1))'];
Ind = [2:N, 1:N, 1:(N-1)];
Jnd = [1:(N-1), 1:N, 2:N];
Mp = sparse(Ind,Jnd,S,N,N);
% negative flux
S = [-fim(2:N)',fim',zeros(1,N-1)];
Ind = [2:N, 1:N, 1:(N-1)];
Jnd = [1:(N-1), 1:N, 2:N];
Mm = sparse(Ind,Jnd,S,N,N);

% compaction right hand side
Shear = CompactionFunction(xgrid,W,phi,Theta,A,nglen,Accumulation,type);

% saturation advection components
SaturationVelocity = ones(N+1,1);
[fadvSp,fadvSm] = AdvectiveFlux(N,SaturationVelocity,U*k(phi));
fdifS = fadvSp-fadvSm;

% water pressure computation
pressure(I) = (-M(I,I))\(dx.*Shear(I)-fdifS(I)+M(I,~I)*pressure(~I));
fpD = -Mp*pressure; fmD = -Mm*pressure;
qp = fpD+fadvSp;
qm = fmD+fadvSm;
end

