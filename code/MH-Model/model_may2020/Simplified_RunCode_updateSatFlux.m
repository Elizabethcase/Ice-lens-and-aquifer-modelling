% script to run code
clearvars; close all; clc;
%%
units = 'mwe';
T = 2; % total simulation time (yr)

%lists options for seasonality & intialization
thetaOpt = {'constant', 'seasonal'}; %(temperature seasonality)
thetaIOpt = {'uniiso', 'tanh', 'gausiso'}; %uniiso = uniform & isothermal; tanh and gausiso both can be customized below
accOpt = {'constant','seasonal'}; %(snow acc. seasonality)
phiOpt = {'exponential', 'gausuni', 'gausexp', 'ice lens uni', 'ice lens exp'}; % porosity initialization
rcOpt = {'constant','seasonal','startstop'};
numRuns = 0;

%sets input for model
ai = 0; %accumulation rate
qb = 0; % surface energy balance
qb_offset = 0; % to allow for avg temp < 0 (or >)
ti = 'tanh'; %temperature initialization uniiso, tanh, gausiso
th = 'constant'; %temperature seasonality - constant or seasonal (highest in summer)
ac = 'constant'; %accumulation seasonality - constant or seasonal (largest in winter)
ph = 'exponential'; %porosity initialization - 'exponential', 'gausuni', 'gausexp', 'ice lens uni', 'ice lens exp'
rv = .1; %[0, 1/40, 10/40] = [0, 1, 10] inch/yr
rt = 'startstop';  %rain seasonality, constant, seasonal, startstop = stops at set point through run
ctype = 'none'; %poreclosure, empirical

write_gif = 1; %0 don't write, 1 write

main(ai, units, qb, qb_offset, T, th, ac, ti, ph, rv, rt, ctype, numRuns, write_gif);

function output = main(ai, units, qb, qb_offset, T, th, ac, ti, ph, rv, rt, ctype, numRuns, write_gif)

global B Stefan U A nglen Pe alpha beta ell Q0 h phi0 rhoi phi_min

% % Physical Parameters % % % % % % % % % % % % % % %
B = 260; % bond number
Stefan = 12; % stefan number
U = 100; % darcy velocity to advection
A = 1; % compaction pressure / viscosity scale
nglen = 1; % Glens law exponent for viscosity
Pe = 11; % Peclet number
alpha = 1; % exponent of capillary pressure
beta = alpha + 1; % saturation flux exponent; default = 2
ell = 20.6; % firn melting lengthscale
Q0 = 200; %energy forcing normalization %W m-2 or kg s-3
h = 14.8; %effective heat transfer coefficient 
rhoi = 917; % ice density
phi_min = 1e-4;
% % % % % % % % % % % % % % % % % % % % % % % % % % %

% % Simulation Parameters % % % % % % % % % % % % % %
plot_amount = 1000; % time between each plot
save_freq = 1000; % frequency at which plots are saved
phi0 = 0.64; % surface porosity
metersofsnow = normalizedaccumulation(ai, units, Q0, phi0);
AccumulationRate = metersofsnow*(1-phi0); % accumulation rate
Qbar = qb; % surface energy flux
CompactionType = ctype;
gifcounter = 1;
% % % % % % % % % % % % % % % % % % % % % % % % % % %

% % Discretization % % % % % % % % % % % % % %
dx = 10^(-2)/2;
dt = 10^(-6);
%dt = 10^(-5);
%% % % % % % % % % % % % % % % % % % % % % % % % % % %


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

% Surface energy balance
if strcmp('seasonal',th)
    %EbarFun = @(tau)Qbar-(Qbar+offset)*cos(2*pi*tau); %Q0 = 200; should this be Qbar-cos(2*pi*tau) ?
    EbarFun = @(tau)Qbar*cos(2*pi*tau)+qb_offset; %Q0 = 200; should this be Qbar-cos(2*pi*tau) ?
elseif strcmp('constant',th)
    EbarFun = @(tau) Qbar;
end

%Temperature of firn column
if strcmp(ti, 'uniiso')
    theta0 = 0;
elseif strcmp(ti, 'tanh')
    start=0;
    range=1; %1 normalized T = 13 K
    ztemp =linspace(pi, -pi, length(xgrid));
    theta0 = range/2*tanh(ztemp) + (start-range/2);
    theta0 = theta0';
elseif strcmp(ti, 'gausiso')
    theta_i = 0;
    aG = -4;
    bG = ell/2;
    cG = 1;
    theta0 = theta_i + aG * exp(-((xgrid .* ell - bG).^2)./(2*cG^2));    
end

% Accumulation
if strcmp('seasonal',ac)
    % Very simple seasonal accumulation (max in winter, min in summer)
    Abar = @(tau)AccumulationRate+AccumulationRate*cos(2*pi*tau) ; 
elseif strcmp('constant',ac)
    Abar = @(tau)AccumulationRate; % Accumulation
end

% Porosity
if strcmp('exponential',ph)
    Lphi = 5; %almost zero by z = 20 (~=ell, quite warm, equiv to juneau, also try increasing)
    phi = phi0 .* exp(-xgrid .* ell/Lphi);
elseif strcmp('gausuni',ph)
    aG = -0.54; %phi (close to ice density)
    bG = ell/2;
    cG = 1;
    phi = phi0 + aG * exp(-((xgrid .* ell - bG).^2)./(2*cG^2));
elseif strcmp('gausexp',ph)
    Lphi = 10; %almost zero by z = 20 (~=ell, quite warm, equiv to juneau, also try increasing)
    phi = phi0 .* exp(-xgrid .* ell/Lphi);
    %currently  just constant density + gaussian, but could also add
    %exponential instead
    aG = -phi(50)*.95; %phi (close to ice density)
    bG = ell/2;
    cG = 1;
    phi = phi + aG * exp(-((xgrid .* ell - bG).^2)./(2*cG^2));
elseif strcmp('uniform',ph)
    phi = phi0 .* ones(N,1); 
elseif strcmp('ice lens uni',ph)
    phi = phi0 .* ones(N,1);
    lensThickness = 2; % # spatial steps;
    zLoc = 5;
    [v,index] = min(abs(zLoc - xcelledges*ell));
    phi(index : index + lensThickness) = 0;
elseif strcmp('ice lens exp',ph)
    Lphi = 10; %almost zero by z = 20 (~=ell, quite warm, equiv to juneau, also try increasing)
    phi = phi0 .* exp(-xgrid .* ell/Lphi);
    lensThickness = 2; % # spatial steps;
    zLoc = 5;
    [v,index] = min(abs(zLoc - xgrid*ell));
    phi(index:index + lensThickness) = 0;
end

%Rainfall at surface
if strcmp(rt, 'constant')
    Rbar = @(tau)rv; % fixed surface water flux (rain)
elseif strcmp(rt,'seasonal')
    Rbar = @(tau)rv-rv*cos(2*pi*tau); %high in summer, zero in winter
elseif strcmp(rt,'startstop')
    Rbar = @(tau)rv; %adjusted in code to go to zero at rstop point
    rstop = 0.5; 
end

% Initial values

pressure0 = -(pc(1)/B)*ones(N,1); % initial pressure
pressure = pressure0;
zs = 0; % zero initial surface height
S0 = zeros(N,1);
%phi = phi0*(1-0.99*exp(-(xgrid-0.5).^2/0.005));
W = 1-phi+phi.*S0; % take in W from above
H =  W.*theta0+Stefan.*phi.*S0; % take in H from above

% Initialize variables
Tsurf = zeros(Nt,1); MeltRate = zeros(Nt,1);
n_plot = 1; RunOff = zeros(Nt,1); RO = zeros(Nt,1);




%%%% loop begins %%%%



for n = 1:Nt
    % Assign values from previous timestep
    W_nm1 = W;
    H_nm1 = H;
    
    if n>1
        Theta_nm1 = Theta;
        S_nm1 = S;
        phi_nm1 = phi;
    else
        [Theta_nm1,phi_nm1,S_nm1] = conversiontotemperature(H_nm1,W_nm1,Stefan);
    end
    
    if abs(n*dt - rstop) < 2.0000e-05
       if strcmp(rt,'startstop')
            Rbar = @(tau) 0;
        end
    end
    
    % Compute Compaction Velocity
    Shear = CompactionFunction(xgrid,W_nm1,phi_nm1,Theta_nm1,A,nglen,Abar(n*dt),CompactionType);
    CompactionVelocity = cumtrapz(xcelledges,[Shear; Shear(end)]);
    SurfaceIceVelocity = (Abar(n*dt)./(1-phi0))- MeltRate(max(n-1,1));
    if SurfaceIceVelocity<0
        SurfaceIceVelocity = (Abar(n*dt)./(1-phi_nm1(1)))- MeltRate(max(n-1,1));
    end
    IceVelocity = SurfaceIceVelocity*ones(N+1,1)-CompactionVelocity;
    
    %% %% %% FLUX CALC STARTS HERE %% %% %%
        
    %% Total Water (Mass) Flux %%
    % Ice Advective flux;
    [FadvIpW,FadvImW] = AdvectiveFlux(N,IceVelocity,W_nm1); %velocity * volume
    
    % Saturation Advective flux - water percolation?
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
    FmW(1) = Abar(n*dt)+Rbar(n*dt); % accumulation and rain :: - RO(max(n-1,1))
    FdifW = FpW-FmW;
    W = W_nm1-(dt/dx)*FdifW;
    
    %% Enthalpy Flux %%
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
    Fm(1)=Stefan*(EbarFun(n*dt)-Theta_nm1(1))+Stefan*Rbar(n*dt); % Enthalpy neumann conditions
    FdifH = Fp-Fm;
    H = H_nm1-(dt/dx)*FdifH;
    
    %%  Saturated water pressure + fluxes %%
    % Compute fully saturated water pressure
    I = S_nm1>=1;
    if I(I)
        [qp,qm,pressure] = FullySaturatedWaterPressure(U,k,pressure,dx,xgrid,N,A,phi_nm1,W_nm1,Theta_nm1,nglen,I,Abar(n*dt),CompactionType);
        if 1
            slocations = find(I(2:(N-1)))+1;
            for i = slocations'            
                if and(S_nm1(i-1)~=1,S_nm1(i)==1) % unsat(i-1)/sat(i)
                    
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
                        
                         % Total Water
                        FpWS(i) = max(qp(i),FpWS(i)); % use maximum
                        FmWS(i+1) = max(qm(i+1),FmWS(i+1)); % use maximum
                        
                        % Enthalpy
                        FpS(i) = max(Stefan*qp(i),FpS(i)); % use maximum
                        FmS(i+1) = max(Stefan*qm(i+1),FmS(i+1)); % use maximum                         
                    end
                        
                elseif and(S_nm1(i)==1,S_nm1(i+1)==1) % sat(i)/sat(i+1)
                    
                    %Total Water
                    FpWS(i) = qp(i); % use q
                    FmWS(i+1) = qm(i+1); % use q
                    
                    % Enthalpy
                    FpS(i) = Stefan*qp(i); % use q
                    FmS(i+1) = Stefan*qm(i+1); % use q
                    
                elseif and(S_nm1(i)==1,S_nm1(i+1)~=1) % sat(i)/unsat(i+1)
  
                    FpWS(i) = max(qp(i),FpWS(i)); % use min or max?
                    FmWS(i+1) = max(qm(i+1),FmWS(i+1)); % use minimum
                        
                    % Enthalpy
                    FpS(i) = max(Stefan*qp(i),FpS(i)); % use minimum
                    FmS(i+1) = max(Stefan*qm(i+1),FmS(i+1)); % use minimum 
                end       
                
            end
            % Fix flux at surface for full saturation
%         elseif 0
%             Ip = or([I(2:N); I(N)],I); Im = or([I(1); I(1:(N-1))],I); % Ipc = [I(2:N); false]; find(Ipc)==Ip;
%             % Total Water
%             FpWS(Ip)=qp(Ip); FmWS(Im)=qm(Im);
%             FpS(Ip)=Stefan*qp(Ip); FmS(Im)=Stefan*qm(Im);
        end
        if I(1) % if saturated at surface
            FmWS(1) = qm(1);
            FmS(1) = Stefan*qm(1);
        elseif I(N) % if saturated at base
            FpWS(N) = qp(N);
            FpS(N) = Stefan*qp(N);
        end
        
        % Total Water
        if I(1) % if saturated at surface
            FpW = FadvIpW+FpWS; FmW = FadvImW+FmWS;
        else
            FpW = FadvIpW + FpWS;
            % Don't doctor FmW(1) if first node not saturated
            % to keep boundary condition imposed earlier
            FmW(2:N) = FadvImW(2:N) + FmWS(2:N);
        end
        % Compute run off
        RO(n) = Abar(n*dt)+Rbar(n*dt)-FmW(1);
        if RO(n)>0
            RunOff(n) = RO(n);
            RunOff_flag=1;
            % FmW(1) = FpW(1);
        else
            FmW(1) = Abar(n*dt)+Rbar(n*dt); % Fixed rain flux
            RunOff_flag = 0;
        end
        
        FdifW = FpW-FmW;
        %W(I) = W_nm1(I)-(dt/dx)*Fdif(I);
        W = W_nm1-(dt/dx)*FdifW;
        
        if RunOff_flag
            W(1)=1;
        end
        
        % Enthalpy
        Fp = FadvIpH+FpS+FpT; Fm = FadvImH+FmS+FmT;
        if ~RunOff_flag
            Fm(1)= Stefan*(EbarFun(n*dt)-Theta_nm1(1)+Rbar(n*dt));
        elseif RunOff_flag
            Fm(1)=Stefan*(EbarFun(n*dt)-Theta_nm1(1)+Rbar(n*dt)-RunOff(n));
        end
        FdifH = Fp-Fm;
        %H(I) = H_nm1(I)-(dt/dx)*Fdif(I);
        H = H_nm1-(dt/dx)*FdifH;
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
        plot(S,xgrid,'r*-','linewidth',2,'DisplayName','Saturation')
        hold on;
        %plot(Theta,xgrid,'y','linewidth',2,'DisplayName','Temperature')
        %plot(H,xgrid,'k','linewidth',2,'DisplayName','Enthalpy')        
        plot(W,xgrid,'b','linewidth',2,'DisplayName','Total Water')
        plot(phi,xgrid,'go-','linewidth',2,'DisplayName','Phi')
        %plot(pressure,xgrid,'m','linewidth',2,'DisplayName','Pressure')
        plot(FpW,xgrid,'c^-','linewidth',1,'DisplayName','FpW')
        plot(FmW,xgrid,'b^-','linewidth',1,'DisplayName','FmW')
        plot(FdifW,xgrid,'^-','Color','#006E77','linewidth',1,'DisplayName','Fdif')
        title(num2str(n*dt))
        set(gca,'fontsize',18,'ydir','reverse')
        axis([-1 2 a b])
        legend();
        drawnow;
        hold off;
        
        %write to gif
        if write_gif==1 
            frame = getframe(1);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            
            if gifcounter == 1
                filename = ['/Users/elizabeth/Documents/projects/Ice-lens-and-aquifer-modelling/code/MH-Model/model_may2020/figures/gif/' datestr(now,'yymmddHHMM') '.gif'];
                imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
                gifcounter = gifcounter + 1;
            else
                imwrite(imind,cm,filename,'gif','WriteMode','append');
            end
        end

       
    end
    

    zs = zs + SurfaceIceVelocity*dt;
    if ~mod(n,save_freq)
        disp(['prog... ' num2str(round(n/save_freq))])
        phiwithz(:,n_plot) = phi;
        Swithz(:,n_plot) = S;
        Thetawithz(:,n_plot)= Theta;
        Wwithz(:,n_plot)= W;
        FdifWwithz(:,n_plot)= FdifW;
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
Wmat = Wwithz;
FdifWmat = FdifWwithz;

tmat = repmat(time,length(xgrid),1);
zmat = ell*repmat(xgrid,1,length(time));

figure('Position',[440 377 560 600])
tiledlayout(5,1)

ax1 = nexttile;
surf(tmat,zmat,phimat,'EdgeColor','none'); view(2); hold on;
hc = colorbar; caxis([0 phi0])
ylabel(hc,'$\phi$','interpreter','latex','fontsize',20)
set(gca,'fontsize',18,'ydir','reverse','layer','top'); grid off;
ylabel('$Z$ (m)','interpreter','latex','fontsize',20)
axis([0 T 0 ell])
plot3(rstop*ones(size(zmat(:,1))), zmat(:,1), ones(size(zmat(:,1)))*phimat(1,1),'-r')
hold off;

ax2 = nexttile;
surf(tmat,zmat,smat,'EdgeColor','none'); view(2);
hc = colorbar; caxis([0 1])
ylabel(hc,'$S$','interpreter','latex','fontsize',20)
set(gca,'fontsize',18,'ydir','reverse','layer','top'); grid off;
ylabel('$Z$ (m)','interpreter','latex','fontsize',20)
axis([0 T 0 ell])

ax3 = nexttile;
deltaT = Q0/h;
surf(tmat,zmat,deltaT*thetamat,'EdgeColor','none'); view(2);
hc = colorbar; caxis auto;
ylabel(hc,'$T$ $^{\circ}$C','interpreter','latex','fontsize',20)
set(gca,'fontsize',18,'ydir','reverse','layer','top'); grid off;
ylabel('$Z$ (m)','interpreter','latex','fontsize',20)
axis([0 T 0 ell])

ax4 = nexttile;
surf(tmat,zmat,FdifWmat,'EdgeColor','none'); view(2);
hc = colorbar; caxis([0 1])
ylabel(hc,'$FdifW$','interpreter','latex','fontsize',20)
set(gca,'fontsize',18,'ydir','reverse','layer','top'); grid off;
ylabel('$Z$ (m)','interpreter','latex','fontsize',20)
axis([0 T 0 ell])

ax5 = nexttile;
plot(tmat,sum(Wmat,1));
ylabel('$\sum W$','interpreter','latex','fontsize',20)
set(gca,'fontsize',18,'layer','top');
xlabel('$t$ (yr)','interpreter','latex','fontsize',20)
axis([0 T min(sum(Wmat,1)) max(sum(Wmat,1))])

sgtitle(sprintf(['T = ' th ', Qbar = ' num2str(qb*Q0) '\nA = ' num2str(ai) ' mwe, ' ac '\n phiMin = ' num2str(min(min(phimat))) '\n RVol = ' num2str(rv)]))

linkaxes([ax1, ax2, ax3, ax4],'xy')
linkaxes([ax1, ax2, ax3, ax4, ax5],'x')

set(gcf,'color','w');


output_new.numRun = numRuns;
output_new.AccumulationRate = ai;
output_new.AccumulationSeasonality = ac;
output_new.AccumulationUnits = units;
output_new.SurfaceEnergyBalance = qb;
output_new.TemperatureInitialShape = th;
output_new.TemperatureSeasonality = th;
output_new.PorosityInitialShape = ph;
output_new.SnowPorosity = phi0;
output_new.RainVolume = rv;
output_new.RainSeasonality = rt;
output_new.CompactionType = ctype;
output_new.dx = dx;
output_new.dt = dt;
output_new.imgFile = [datestr(now,'yymmddHHMM') '.fig'];
output_new.date = datetime('now','Format','yyMMdd');

answer = questdlg('Would you to save the figure?', ...
	'SaveFig', ...
	'Yes','No',...
    'Yes');
% Handle response
switch answer
    case 'Yes'
        %save figures
        savefig(['/Users/elizabeth/Documents/projects/Ice-lens-and-aquifer-modelling/code/MH-Model/model_may2020/figures/fig/' datestr(now,'yymmddHHMM')]);
        saveas(gcf, ['/Users/elizabeth/Documents/projects/Ice-lens-and-aquifer-modelling/code/MH-Model/model_may2020/figures/png/' datestr(now,'yymmddHHMM') '.png']);
        %concatenate & save output as table
        load('Output_MHMay.mat','output');
        output2 = struct2table(output_new);
        output_join = [output;output2];
        clear output
        output = output_join;
        save('/Users/elizabeth/Documents/projects/Ice-lens-and-aquifer-modelling/code/MH-Model/model_may2020/Output_MHMay.mat','output');
    case 'No'
end

%%% create fresh output table
% varNames=["numRun", "AccumulationRate", "AccumulationSeasonality", "AccumulationUnits", "SurfaceEnergyBalance",..."TemperatureInitialShape", "TemperatureSeasonality","PorisityInitialShape", "SnowPorosity","RainVolume", "RainSeasonality","CompactionType","dx", "dt","imgFile", "date"];
% sz = [1 length(varNames)]
% varTypes=["double","double","string","string","double","string","string","string","double","double","string","string","double","double","string","datetime"]
% output = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames)

end