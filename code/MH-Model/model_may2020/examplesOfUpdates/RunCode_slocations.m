% script to run code
clearvars; close all; clc;

countp = [];
countm = [];
% % Simulation Parameters % % % % % % % % % % % % % %
T = 1; % total simulation time (yr)
plot_amount = 1000; % time between each plot
save_freq = 100; % frequency at which plots are saved
phi0 = 0.64; % surface porosity
AccumulationRate = 0.0*(1-phi0); % accumulation rate
Qbar = -0.8; % surface energy flux
type = 'none';
% % % % % % % % % % % % % % % % % % % % % % % % % % %

% % Discretization % % % % % % % % % % % % % % % % % 
dx = 10^(-2); dt = 10^(-5);
% % % % % % % % % % % % % % % % % % % % % % % % % % %

% % Physical Parameters % % % % % % % % % % % % % % %
B = 260; % bond number
Stefan = 12; % stefan number = Latent heat/ (cp * dT) where Lh = 334,000 m2s-2; cp = 2050 m2s-2K-1; dT = 13.5K = Q0/h
U = 100; % darcy velocity to advection
A = 1; % compaction pressure / viscosity scale
nglen = 1; % Glens law exponent for viscosity
Pe = 11; % Peclet number
alpha = 1; % exponent of capillary pressure
beta = 2; % saturation flux exponent
ell = 20.6; % firn melting lengthscale
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
EbarFun = @(tau)0*(Qbar+cos(2*pi*tau));

% Initial values
Rbar = 0.1; % fixed surface water flux (rain)
pressure0 = -(pc(1)/B)*ones(N,1); % initial pressure
zs = 0; % zero initial surface height
S0 = zeros(N,1);
theta0=0;
phi = phi0 * (exp(-(xgrid).^2/.3)); %0.8*phi0*(1 - 0.8*xgrid); 
W = 1-phi+phi.*S0; % take in W from above
H =  W.*theta0+Stefan.*phi.*S0; % take in H from above

% Initialize variables
Tsurf = zeros(Nt,1); MeltRate = zeros(Nt,1);
n_plot = 1; RunOff = zeros(Nt,1); RO = zeros(Nt,1);
for n = 1:Nt
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
    
    % Total water saturation flux (liquid advective + diffusive)
    FpWS = FadvSp+FpD;
    FmWS = FadvSm+FmD;
    
    % Total fluxes (liquid + carried by ice)
    FpW = FadvIpW+FpWS; FmW = FadvImW+FmWS;
    FmW(1) = Abar(n*dt)+Rbar; % accumulation and rain :: - RO(max(n-1,1))
    Fdif = FpW-FmW;
    W = W_nm1-(dt/dx)*Fdif;
    
    % % Enthalpy % %
    % Ice Advective flux;
    [FadvIpH,FadvImH] = AdvectiveFlux(N,IceVelocity,H_nm1);
    
    % Saturation Advective flux;
    SaturationVelocity = ones(N+1,1);
    [fadvSp,fadvSm] = AdvectiveFlux(N,SaturationVelocity,k(phi_nm1).*kr(S_nm1)); %same as above?
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
        if 1
            slocations = find(I(2:(N-1)))+1; % ignores saturation at surface or base (may not need to ignore saturation at base?); 
            for i = slocations % *facepalm*               
                if and(S_nm1(i-1)~=1,S_nm1(i)==1) % unsat(i-1)/sat(i)
                    
                    % Total Water
                    
                    if min(qp(i-1),FpWS(i-1)) == qp(i-1)
                        disp('min - ps')
                        disp(['qp = ' num2str(qp(i-1)) ' | up = ' num2str(FpWS(i-1))]);
                    else 
                        disp('min - pu')
                    end
                    if min(qm(i),FmWS(i)) == qm(i)
                        disp('min - ms')
                        disp(['qm = ' num2str(qm(i)) ' | um = ' num2str(FmWS(i))]);                        
                    else 
                        disp('min - mu')
                    end
                    
                    FpWS(i-1) = min(qp(i-1),FpWS(i-1));
                    FmWS(i) = min(qm(i),FmWS(i));  
                    
                    % Enthalpy
                    FpS(i-1) = min(Stefan*qp(i-1),FpS(i-1)); % use minimum
                    FmS(i) = min(Stefan*qm(i),FmS(i)); % use minimum
                    
                    if and(S_nm1(i)==1,S_nm1(i+1)==1) % sat(i)/sat(i+1)
                       
                        % Total Water
                        FpWS(i) = qp(i); % use q
                        FmWS(i+1) = qm(i+1); % use q
                        % Enthalpy
                        FpS(i) = Stefan*qp(i); % use q
                        FmS(i+1) = Stefan*qm(i+1); % use q    
                        
                    elseif and(S_nm1(i)==1,S_nm1(i+1)~=1) % sat(i)/unsat(i+1)

                        if max(qp(i),FpWS(i)) == qp(i)
                            disp('max - ps')
                            disp(['qp = ' num2str(qp(i)) ' | up = ' num2str(FpWS(i))]);
                        else 
                            disp('max - pu')
                        end

                        if max(qm(i+1),FmWS(i+1)) == qm(i+1)
                            disp('max - ms')
                            disp(['qm = ' num2str(qm(i+1)) ' | um = ' num2str(FmWS(i+1))]);   
                        else 
                            disp('max - mu')     
                        end
                        disp('--')
                         % Total Water
                        FpWS(i) = max(qp(i),FpWS(i)); % use maximum
                        FmWS(i+1) = max(qm(i+1),FmWS(i+1)); % use maximum
                        % Enthalpy
                        FpS(i) = max(Stefan*qp(i),FpS(i)); % use maximum
                        FmS(i+1) = max(Stefan*qm(i+1),FmS(i+1)); % use maximum                         
                    end
                        
                elseif and(S_nm1(i)==1,S_nm1(i+1)==1) % sat(i)/sat(i+1)
                    FpWS(i) = qp(i); % use q
                    FmWS(i+1) = qm(i+1); % use q
                    % Enthalpy
                    FpS(i) = Stefan*qp(i); % use q
                    FmS(i+1) = Stefan*qm(i+1); % use q
                    
                elseif and(S_nm1(i)==1,S_nm1(i+1)~=1) % sat(i)/unsat(i+1)
                    % Total Water

                        if max(qp(i),FpWS(i)) == qp(i)
                            disp('max - ps')
                            disp(['qp = ' num2str(qp(i)) ' | up = ' num2str(FpWS(i))]);
                        else 
                            disp('max - pu')
                        end

                        if max(qm(i+1),FmWS(i+1)) == qm(i+1)
                            disp('max - ms')
                            disp(['qm = ' num2str(qm(i+1)) ' | um = ' num2str(FmWS(i+1))]);   
                        else 
                            disp('max - mu')     
                        end   
                    FpWS(i) = min(qp(i),FpWS(i)); % use maximum
                    FmWS(i+1) = min(qm(i+1),FmWS(i+1)); % use maximum
                    % Enthalpy
                    FpS(i) = min(Stefan*qp(i),FpS(i)); % use maximum
                    FmS(i+1) = min(Stefan*qm(i+1),FmS(i+1)); % use maximum 
                end       
                
            end
            % Fix flux at surface for full saturation
        elseif 0
            Ip = or([I(2:N); I(N)],I); Im = or([I(1); I(1:(N-1))],I); % Ipc = [I(2:N); false]; find(Ipc)==Ip;
            % Total Water
            FpWS(Ip)=qp(Ip); FmWS(Im)=qm(Im);
            FpS(Ip)=Stefan*qp(Ip); FmS(Im)=Stefan*qm(Im);
        end
        if I(1)
            FmWS(1) = qm(1);
            FmS(1) = Stefan*qm(1);            
        elseif I(N)
            FpWS(N) = qp(N);
            FpS(N) = Stefan*qp(N);            
        end
%         if I(1) ~= 1
%             FmWS(1) = FmWS(1);
%             FmS(1) = FmS(1);
%         elseif I(1) == 1
%             FmWS(1) = qm(1);
%             FmS(1) = Stefan*qm(1);            
%         elseif I(N) ~= 1
%             FpWS(N) = FpWS(N);
%             FpS(N) = FpS(N);
%         elseif I(N) == 1
%             FpWS(N) = qp(N);
%             FpS(N) = Stefan*qp(N);            
%         end
        
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
       %[Theta,phi,S] = conversiontotemperature(H,W,Stefan); %%%%% turn off for debug

    else
        pressure = pressure0;
    end
    disp('***')
    %%%% end of (full) saturation %%%%
    
    [Theta,phi,S] = conversiontotemperature(H,W,Stefan); %%%%% turn off for debug
    [H,W] = conversiontoenthalpy(Theta,phi,S,Stefan); %%%% turn off for debug
    
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
        plot(H,xgrid,'k','linewidth',2)
        hold on;
        plot(Theta,xgrid,'y','linewidth',2)
        plot(S,xgrid,'r','linewidth',2)
        plot(W,xgrid,'b','linewidth',2)
        plot(phi,xgrid,'g','linewidth',2)
        plot(pressure,xgrid,'m','linewidth',2)
        plot(FmWS,xgrid,'c','linewidth',2)
        plot(FpWS,xgrid,'c--','linewidth',2)
        title(num2str(n*dt))
        set(gca,'fontsize',18,'ydir','reverse')
        axis([-1 2 a b])
        drawnow;
        hold off;
    end
    
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

figure(2)
subplot(3,1,1)
surf(tmat,zmat,phimat,'EdgeColor','none'); view(2);
hc = colorbar; caxis([0 phi0])
ylabel(hc,'$\phi$','interpreter','latex','fontsize',20)
set(gca,'fontsize',18,'ydir','reverse','layer','top'); grid off;
ylabel('$Z$ (m)','interpreter','latex','fontsize',20)
axis([0 T 0 ell])

subplot(3,1,2)
surf(tmat,zmat,smat,'EdgeColor','none'); view(2);
hc = colorbar; caxis([0 1])
ylabel(hc,'$S$','interpreter','latex','fontsize',20)
set(gca,'fontsize',18,'ydir','reverse','layer','top'); grid off;
ylabel('$Z$ (m)','interpreter','latex','fontsize',20)
axis([0 T 0 ell])

subplot(3,1,3)
deltaT = 200/14.8;
surf(tmat,zmat,deltaT*thetamat,'EdgeColor','none'); view(2);
hc = colorbar; caxis([-2,0])
ylabel(hc,'$T$ $^{\circ}$C','interpreter','latex','fontsize',20)
set(gca,'fontsize',18,'ydir','reverse','layer','top'); grid off;
xlabel('$t$ (yr)','interpreter','latex','fontsize',20)
ylabel('$Z$ (m)','interpreter','latex','fontsize',20)
axis([0 T 0 ell])
