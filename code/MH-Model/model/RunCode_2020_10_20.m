% script to run code
%clearvars; close all; clc;

% AArray = 1.2; %linspace(4,10,6); units = 'mwe';
% QArray = -0.5; %-linspace(0.5, 1, 5);
% 
% %A,Q = {-0.5,1}, {-.5,2}, {-.667, 1} creates ice lenses
% 
Q0 = 200; %energy forcing %W m-2 or kg s-3
t0 = 3.15*10^7; % seconds per year s
rhoI = 917; % ice denisty kg m-3
L = 334000; %latent heat of water; m2 s-2
phi0 = 0.64;
% 
% 
NumRuns = 0;
% 
% for iA = 1:length(AArray)
%     AIn = normalizedaccumulation(AArray(iA), units, Q0, t0, rhoI, L, phi0); 
%     for iQ = 1:length(QArray)
%         NumRuns = NumRuns+1;
%         disp(['Qbar: ' num2str(QArray(iQ)) ' W m-2;  A: ' num2str(AArray(iA)) ' m.w.e.'])
%         main(AIn, QArray(iQ), NumRuns)
%     end
% end

QbarIn = -0.5;
mwe = 4;
units='mwe';
AIn = normalizedaccumulation(mwe, units, Q0, t0, rhoI, L, phi0);
main(AIn, QbarIn, NumRuns)

function main(AIn, QbarIn, NumRuns)

% % Simulation Parameters % % % % % % % % % % % % % %
T = 10; % total simulation time (yr)
plot_amount = 10000; % time between each plot s
save_freq = 100; % frequency at which plots are saved
phi0 = 0.64; % surface porosity
metersofsnow = AIn;
AccumulationRate = metersofsnow*(1-phi0); %metersofsnow*(1-phi0); % m snow / yr 
Qbar = QbarIn; %; % surface energy flux normalized by ? 
type = 'empirical'; %none, poreclosure
% % % % % % % % % % % % % % % % % % % % % % % % % % %

% % Discretization % % % % % % % % % % % % % % % % % %
dx = 10^(-2); 
dt = 10^(-4);   
% % % % % % % % % % % % % % % % % % % % % % % % % % %

% % Physical Parameters % % % % % % % % % % % % % % %
B = 260; % bond number
Stefan = 12; % stefan number
U = 100; % darcy velocity to advection
A = 1; % compaction pressure / viscosity scale
nglen = 1; % Glens law exponent for viscosity
Pe = 11; % Peclet number
alpha = 1; % exponent of capillary pressure
beta = 2; % saturation flux exponent
ell = 20.6; % firn melting lengthscale
t0 = 3.15*10^7; % seconds per year s
% % % % % % % % % % % % % % % % % % % % % % % % % % %

% Mesh Information
% domain from from x=a to x=b
a = 0;
b = 1;
N = round((b-a)/dx); % Number of grid cells
xcelledges = linspace(a,b,N+1)'; % Cell edges
xgrid = (xcelledges(1:N)+xcelledges(2:(N+1)))/2; % Cell centers

% Time step information
Nt = round(T/dt); % Number of timesteps

% Functions
k = @(p) p.^(3); % Simple Carmen-Kozeny
kr = @(s) s.^beta; % relative permeability with saturation
pc = @(s) s.^(-alpha); % capillary pressure
kr_pc_prime = @(s) -alpha.*(s.^(beta-(alpha+1))); % combined function
Abar = @(tau)AccumulationRate; % Accumulation
EbarFun = @(tau)Qbar-cos(2*pi*tau); %Q0 = 200

% Initial values
Rbar = 8.05 * 10^-10; %0; % fixed surface water flux (rain)
pressure0 = -(pc(1)/B)*ones(N,1); % initial pressure
zs = 0; % zero initial surface height
W = (1-phi0)*ones(N,1); % take in W from above
H =  W.*Qbar; % take in H from above

% Initialize variables
Tsurf = zeros(Nt,1); MeltRate = zeros(Nt,1);
n_plot = 1; RunOff = zeros(Nt,1); RO = zeros(Nt,1);

holdOn = 0; %allows warming to stick around for longer
for n = 1:Nt
    
    % Ebarfun has warming spike in winter
    % dt = 10^-4 ~ 2.89 days.. would be a long warming spike.. may need to
    % decrease dt
    holdOn = holdOn + 1;
    if holdOn > 20
    EbarFun = @(tau)Qbar-cos(2*pi*tau);
        if EbarFun(n*dt) < Qbar % if we are in coldest half of year (E < average)
            disp('winter')
            if rand(1) < .02 % warming spikes 2 percent of the time
                disp('triggered')
                EbarFun = @(tau) Qbar+1; %
                holdOn = 0;
            end
        end
    end

    
    % Assign values from previous timestep
    W_nm1 = W;
    H_nm1 = H;
    % convert to Theta, phi, and S
    [Theta_nm1,phi_nm1,S_nm1] = conversiontotemperature(H_nm1,W_nm1,Stefan);
    
    % % Theta is temperature % phi is porosity % S is saturation % %
    
    % Compute Compaction Velocity
    Shear = CompactionFunction(xgrid,W_nm1,phi_nm1,Theta_nm1,A,nglen,Abar(n*dt),type);
    CompactionVelocity = cumtrapz(xcelledges,[Shear; Shear(end)]);
    SurfaceIceVelocity = (Abar(n*dt)./(1-phi0))- MeltRate(max(n-1,1));
    if SurfaceIceVelocity<0
        SurfaceIceVelocity = (Abar(n*dt)./(1-phi_nm1(1)))- MeltRate(max(n-1,1));
    end
    IceVelocity = SurfaceIceVelocity*ones(N+1,1)-CompactionVelocity;
    
    % % Total Water % %
    % Ice Advective flux;
    [FadvIpW,FadvImW] = AdvectiveFlux(N,IceVelocity,W_nm1);
    
    % Saturation Advective flux;
    SaturationVelocity = ones(N+1,1);
    [fadvSp,fadvSm] = AdvectiveFlux(N,SaturationVelocity,k(phi_nm1).*kr(S_nm1));
    FadvSp = U*fadvSp;
    FadvSm = U*fadvSm;
    
    % Diffusive flux
    [FpD,FmD] = SaturationDiffusiveFlux(@(x,y)k(x).*kr_pc_prime(y),phi_nm1,S_nm1,dt,dx,U/B);
    
    % No diffusive flux boundary condition
    FpD(end)=0;
    
    % Total water saturation flux
    FpWS = FadvSp+FpD;
    FmWS = FadvSm+FmD;
    
    % Total fluxes
    FpW = FadvIpW+FpWS; FmW = FadvImW+FmWS;
    FmW(1) = Abar(n*dt)+Rbar; % accumulation and rain :: - RO(max(n-1,1))
    Fdif = FpW-FmW;
    W = W_nm1-(dt/dx)*Fdif;
    
    % % Enthalpy % %
    % Ice Advective flux;
    [FadvIpH,FadvImH] = AdvectiveFlux(N,IceVelocity,H_nm1);
    
    % Saturation Advective flux;
    SaturationVelocity = ones(N+1,1);
    [fadvSp,fadvSm] = AdvectiveFlux(N,SaturationVelocity,k(phi_nm1).*kr(S_nm1));
    FadvSp = U*Stefan*fadvSp;
    FadvSm = U*Stefan*fadvSm;
    
    % Enthalpy/Saturation Diffusive flux
    [FpE,FmE] = EnthalpyDiffusiveFlux(@(x,y)k(x).*kr_pc_prime(y),Theta_nm1,phi_nm1,S_nm1,dt,dx,U/B,Stefan);
    
    % Temperature Diffusive flux
    [FpT,FmT] = TemperatureDiffusiveFlux(W_nm1,Theta_nm1,dt,dx,1/Pe);
    
    % No diffusive flux boundary condition
    FpT(end) = 0; FpE(end) = 0;
    
    % Total saturation flux
    FpS = FpE+FadvSp; FmS = FmE+FadvSm;
    
    % Total Saturation and Temperature fluxes
    Fp = FadvIpH+FpS+FpT; Fm = FadvImH+FmS+FmT;
    Fm(1)=Stefan*(EbarFun(n*dt)-Theta_nm1(1))+Stefan*Rbar; % Enthalpy neumann conditions
    Fdif = Fp-Fm;
    H = H_nm1-(dt/dx)*Fdif;
    
    % Compute fully saturated water pressure
    I = S_nm1>=1;
    if I(I)
        [qp,qm,pressure] = FullySaturatedWaterPressure(U,k,pressure0,dx,xgrid,N,A,phi_nm1,W_nm1,Theta_nm1,nglen,I,Abar(n*dt),type);
        slocations = find(I(2:(N-1)))+1;
        for i = slocations
            if and(S_nm1(i-1)~=1,S_nm1(i)==1) % unsat(left)/sat(right)
                % Total Water
                FpWS(i-1) = min(qp(i-1),FpWS(i-1)); % use minimum
                FmWS(i) = min(qp(i-1),FpWS(i-1)); % use minimum
                % Enthalpy
                FpS(i-1) = min(Stefan*qp(i-1),FpS(i-1)); % use minimum
                FmS(i) = min(Stefan*qp(i-1),FpS(i-1)); % use minimum
            elseif and(S_nm1(i)==1,S_nm1(i+1)==1) % sat(left)/sat(right)
                % Total Water
                FpWS(i) = qp(i); % use q
                FmWS(i+1) = qp(i); % use q
                % Enthalpy
                FpS(i) = Stefan*qp(i); % use q
                FmS(i+1) = Stefan*qp(i); % use q
            elseif and(S_nm1(i)==1,S_nm1(i+1)~=1) % sat(left)/unsat(right)
                % Total Water
                FpWS(i) = max(qp(i),FpWS(i)); % use maximum
                FmWS(i+1) = max(qp(i),FpWS(i)); % use maximum
                % Enthalpy
                FpS(i) = max(Stefan*qp(i),FpS(i)); % use maximum
                FmS(i+1) = max(Stefan*qp(i),FpS(i)); % use maximum
            end
        end
        % Fix flux at surface for full saturation
        if I(1)
            FmWS(1) = qm(1);
            FmS(1) = Stefan*qm(1);
        elseif I(N)
            FpWS(N) = qp(N);
            FpS(N) = Stefan*qp(N);
        end
        
        % Total Water
        if I(1)
            FpW = FadvIpW+FpWS; FmW = FadvImW+FmWS;
        else
            FpW = FadvIpW + FpWS;
            % Don't doctor FmW(1) if first node not saturated
            % to keep boundary condition imposed earlier
            FmW(2:N) = FadvImW(2:N) + FmWS(2:N);
        end
        % Compute run off
        RO(n) = Abar(n*dt)+Rbar-FmW(1);
        if RO(n)>0
            RunOff(n) = RO(n);
            RunOff_flag=1;
            % FmW(1) = FpW(1);
        else
            FmW(1) = Abar(n*dt)+Rbar; % Fixed rain flux
            RunOff_flag = 0;
        end
        Fdif = FpW-FmW;
        W = W_nm1-(dt/dx)*Fdif;
        if RunOff_flag
            W(1)=1;
        end
        
        % Enthalpy
        Fp = FadvIpH+FpS+FpT; Fm = FadvImH+FmS+FmT;
        if ~RunOff_flag
            Fm(1)= Stefan*(EbarFun(n*dt)-Theta_nm1(1)+Rbar);
        elseif RunOff_flag
            Fm(1)=Stefan*(EbarFun(n*dt)-Theta_nm1(1)+Rbar-RunOff(n));
        end
        Fdif = Fp-Fm;
        H = H_nm1-(dt/dx)*Fdif;
    else
        pressure = pressure0;
    end
    [Theta,phi,S] = conversiontotemperature(H,W,Stefan);
    [H,W] = conversiontoenthalpy(Theta,phi,S,Stefan);
    
    % Compute surface melt rate
    Tsurf(n) = Theta(1);
    DiffusiveFlux = -(2/(Pe*dx))*(((1./W(1))+(1./W(2))).^(-1))*(Theta(2)-Theta(1));
    MR = (EbarFun(n*dt)-(DiffusiveFlux/Stefan))./(1-phi(1));
    ThetaSurface = (3/2)*Theta(1) - (1/2)*Theta(2);
    %     if and(Theta(1)==0,MR>0)
    if and(ThetaSurface>=0,MR>0)
        MeltRate(n) = MR;
    else
        MeltRate(n)=0;
    end
    
    if ~mod(n,plot_amount)
        plot(H,xgrid,'k','linewidth',2, 'DisplayName', 'Enthalpy')
        hold on;
        plot(Theta,xgrid,'y','linewidth',2,'DisplayName','Temperature')
        plot(S,xgrid,'r','linewidth',2,'DisplayName','Saturation')
        plot(W,xgrid,'b','linewidth',2,'DisplayName','Total Water')
        plot(phi,xgrid,'g','linewidth',2,'DisplayName','porosity')
        plot(pressure,xgrid,'m','linewidth',2,'DisplayName','water pressure')
        plot(FmWS,xgrid,'c','linewidth',2,'DisplayName','water flux')
        title(num2str(n*dt))
        set(gca,'fontsize',18,'ydir','reverse')
        axis([-1 2 a b])
        xlabel('Normalized parameter value')
        ylabel('Depth')
        legend();
        drawnow;
        hold off;
        
    end
%     
    zs = zs + SurfaceIceVelocity*dt;
    if ~mod(n,save_freq)
        phiwithz(:,n_plot) = phi;
        Swithz(:,n_plot) = S;
        Thetawithz(:,n_plot)= Theta;
        MeltRateSave(n_plot) = MeltRate(n);
        RunOffSave(n_plot) = RunOff(n);
        time(n_plot) = n*dt;
        SIV(n_plot) = IceVelocity(1); % surface ice velocity
        BIV(n_plot) = IceVelocity(end); % bottom ice velocity
        BottomIceFlux(n_plot) = FadvImW(end); 
        IceFlux(n_plot) = FadvImW(1);
        WaterFlux(n_plot) = FmWS(1);
        BottomWaterFlux(n_plot) = FpWS(end);
        TotalFlux(n_plot) = FmW(1);
        n_plot = n_plot+1;
    end
    LiquidWater(n) = trapz(xgrid,phi.*S);
    TotalIce(n) = trapz(xgrid,1-phi);
end

% Save solution
tvec = time;
MR = MeltRateSave;
RO = RunOffSave;
phimat = phiwithz;
smat = Swithz;
thetamat = Thetawithz;
drainage = -BIV.*phimat(end,:).*smat(end,:) + BottomWaterFlux;

tmat = repmat(time,length(xgrid),1);
zmat = ell*repmat(xgrid,1,length(time));

%save([num2str(NumRuns) '_type_' type '_' datestr(now,'mmddyyHH')])
%save(['Q_' num2str(Qbar), '_A_' num2str(A) '_type_' type '_' datestr(now,'mmddyyHH')])
zmax = ell;

figure()
subplot(3,1,1)
surf(tmat,zmat,phimat,'EdgeColor','none'); view(2);
hc = colorbar; %set(hc,'ylim',[0, phi0]);%[min(min(phimat)) phi0])
caxis([0,phi0])
ylabel(hc,'$\phi$','interpreter','latex','fontsize',20)
set(gca,'fontsize',18,'ydir','reverse','layer','top'); grid off;
ylabel('$Z$ (m)','interpreter','latex','fontsize',20)
axis([0 T 0 zmax])

subplot(3,1,2)
surf(tmat,zmat,smat,'EdgeColor','none'); view(2);

% if max(max(smat)) > 0
%     hc = colorbar; set(hc,'ylim',[0 max(max(smat))])
% else
hc = colorbar; 
%set(hc,'ylim',[0 1]);
caxis([0 1]);
%end
ylabel(hc,'$S$','interpreter','latex','fontsize',20)
set(gca,'fontsize',18,'ydir','reverse','layer','top'); grid off;
ylabel('$Z$ (m)','interpreter','latex','fontsize',20)
axis([0 T 0 zmax])

subplot(3,1,3)
deltaT = 200/14.8;
surf(tmat,zmat,deltaT*thetamat,'EdgeColor','none'); view(2);
hc = colorbar; 
%set(hc,'ylim',[-20,0])
caxis([-20,0]);
ylabel(hc,'$T$ $^{\circ}$C','interpreter','latex','fontsize',20)
set(gca,'fontsize',18,'ydir','reverse','layer','top'); grid off;
xlabel('$t$ (yr)','interpreter','latex','fontsize',20)
ylabel('$Z$ (m)','interpreter','latex','fontsize',20)
axis([0 T 0 zmax])

sgtitle(['Qbar = ' num2str(Qbar) ', A = ' num2str(metersofsnow) ', type = ' type])

savefig([num2str(NumRuns) '_type_' type '_' datestr(now,'mmddyyHH')]);

end