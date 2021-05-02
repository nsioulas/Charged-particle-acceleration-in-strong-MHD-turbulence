%   v. 170216_0040
% Monte Carlo model for Fermi acceleration
%
%   The setup includes three model types
%       * imodel 0: no scattering
%       * imodel 1: original Fermi-like model
%       * imodel 2: Current sheets \Delta W  = |q| L E
%       * imodel 3: Current sheets (W,v)-dependent \Delta W
%       * imodel 12: combination of 1 and 2
%       * imodel 13: combination of 1 and 3
%   with three types of boundary conditions
%       * iBCs 0: open
%       * iBCs 1: periodic
%       * iBCs 2: closed with reinserting
%   in combination with collisions
%
% Adding timing
% Monitoring kicks
% ------------------------------------------------------------------------
clear all
rng('default')
%% User defined parameters
% Model
imodel = 2;            % model type []
                        %   0: no scattering
                        %   1: Fermi-like
                        %   2: Current sheets \Delta W  = |q| L E
                        %   3: Current sheets, (W,v)-dependent with (V/c)^2 term
                        %  12: combination of 1 and 2
                        %  13: combination of 1 and 3

PdW = 0;              %      Probability of model 1

iBCs = 0;               % boundary conditions type []
                        %   0: open
                        %   1: periodic
                        %   2: closed with reinserting
ndims = 3;              % number of dimensions

iCollFlag = 0;          % collisions model
                        %   0: No collisions
                        %   1: Lenard-Bernstein model

ldWTerm_cosPhi  = true; % dW = [...] * ([...] + V*v * cos(phi))
                        %                       maybe combined with next
ldWTerm_dBB     = false;% dW = [...] * ([...] + V*v * dB / B)
                        %                       maybe combined with previous
idWTerm_W       = true; % 1: dW = 2*GammaV * WiTot * ([...])
                        % 0: dW = 2*GammaV * Wrest * ([...])
idWTerm_const   = 1;    % 1: dW = [...] * (V^2/c^2 + [...])
%%  Particles                       
nP   = 1e6;             % Number of particles
N1   = 2*1e3;           % Number of iterations for each particle   (You have to figure out the number you need in order...) 
kT0  = 100;             % initial energy [eV]
lIon = true;            % type of particles: ions (T) or electrons (F)
Zion = 1;               % particles' atomic number (used for ions only) []
v0thresFactor = 0;      % inital speed threshold; i.e. allow only particles
                        %   with speed above vth * v0thresFactor
WmonitFactor = 0;       % monitor energy threshold; i.e. monitor only particles
                        %   with energy above WmonitFactor * Wth
%% Monitoring
% Set either ntM or tfin to 0 for no monitoring
ltfinStop    = true;    % stop when reaching t_final; usefull when iBCs == 0
tinit        = 0;       % initial time [sec]
tfin         = 2;       % final monitoring time [sec]
t_init_monit = 1.e-6;
t_final=1.e2;           %final simulation time [sec]
ntM = 3;                % number of monitors []
% Full tracking
%lTScomplete = true;    % Keep tracking of all scattering related values:  %   t, W, dW, dtScatt, kicks, and model (when random)

%% INITIAL POSITION AND ENERGY OF PARTICLES
 l0centre = false;       % start from centre (T) or from random point (F)

 lW0const = false;       % constant (T) initial energy or from distribution (F)       

%% Enable NON-RECONNECTING UCS
  
  Non_reconn_UCS=0;       % If set to 1 then you have both reconnecting and non reconn UCS         NS
                          % with probability (1) for P(rec)=P_reconnect
                          % (2)for non-reconnecting  P(non- rec)= 1-P_reconnect
                       
  
    P_reconnect=1;        % Probability to encounter a reconnecting current sheet   %NS
                          % P=1: only reconnecting
                          % P=0: non-reconnecting 
                         
  if ( Non_reconn_UCS) &&  (imodel==12 ) 
      
      Non_reconn_UCS=0;  % We dont want to have both non-rec and Alfvenic scatterers   NS
      P_reconnect=1;
  end
%% DECAYING MAGNETIC FIELD

 magnetic_field_decay = false;   % if true we have a gradual conversion of UCSs to AS
                                 % with  a rate dictated by the desired
                                 % function                                   
PdW_init = 0.095;                % initial probability for model 1

t_decay  = 1e-11;                % decaying magnetic field                                                         
 if (imodel==1) || (imodel==2 && (P_reconnect==1 || Non_reconn_UCS==0))
     
    magnetic_field_decay=false;      
 end 
 %% TURN BOX TO RECTANGLE                             
square_to_rectangle = 0;    % If set to 1 then the box becomes a rectangle of dimensions            NS
                            % dx=dy=L and dz=Lz


L = 1.e9;                   % length of the box size [cm]
Lz = 1.e5;                  % If square_to_rectangle=0 then Lz's value does
                            % not matter                      
 %% AT T=0 all paticles on the fractal or not?                      
                       
 energ_t0 = 0;             %if 1 :  energization event for each particle at t=0sec      NS
                           % if 0 : Particle firt performs rand step then 
                           % energization event                         
   
 Dw_at_t0_is_0_or_Dw = 1;  % because the particle has to begin at the fractal           NS
                           % but we may not want to energize it at t=0  
                           % so if 0=> no energization at t=0, altough on fractal                          
%% PROBABLITY DENSITY FUNCTION P(Dr)  -> MONTE CARLO                         
upper_step_limit = true;   % If (F) dr can theoretically take values up to infinity        NS
                           %  i.e. escape from box without second scattering
 P_law_index     = 1.2;    % P(dr)~dr^(P_law_index)
 
 Dr_min          = 10^2;   % Define minimum step size [cm]
 Dr_max          =L;       % Define maximum step size [cm]   
 
%% Scatterers
B       = 100;             % magnetic field [Gauss]
Vfactor = 1;               % Cloud speed = (Alfven velocity) x Vfactor []
n0      = 1.e9;            % number density [cm^(-3)]

%% Box size
L2   = L/2;
dmin = -L2;    % To check if the particle has left the box
dmax =  L2;

dmin_z = -L2;  % to make the code more simple (for the square_to_rectangle case) we set dmin_z
dmax_z =L2;    % seperatly alhtough it may take tha same values as dx,dy(in the cube case)

if square_to_rectangle
  
   L2z    = Lz/2;
   dmin_z = -L2z;
   dmax_z = L2z;
end
%% Model specific parameters
%% For reconnecting UCS
BPLa = 5/3;         % PL spectrum: p(??)~??^(BPLa)
BLim1 = 1.e-5;      % Lower limit
BLim2 = 1.e2;       % Upper limit
% UCS; model 2
idtCS = 0;          % 0: consider instantaneous scattering; 1: not
%   Beff
lBdepend = true;    % Eeff & Leff depend on Beff

%   Leff
LeffLim1 = 1.e2;    % Lower limit (or constant value, when Leff is a free parameter)
LeffLim2 = 1.e5;    % Upper limit
lLeffrand = false;  % When Leff is a free parameter: Uniform distr. (T) or constant (F)

%% For NON reconnecting UCS
second_way=false;       %False: Just ignore collision
%%   Leff
LeffLim1_non = 1.e2;    % Lower limit (or constant value, when Leff is a free parameter)
LeffLim2_non = 1.e5;    % Upper limit
lLeffrand_non = false;  % When Leff is a free parameter: Uniform distr. (T) or constant (F)
%%  Eeff
T0cs = 1.e6;            % Temperature (used for Dreicer field) [K]
%Power-Law distribution
lEplaw_non =true;       % Power-law distribution (T) or constant Eeff for imodel = 2
EPLa_non = 5/3;         % PL spectrum: x^{EPLa}
% The following limits will be relative to the Dreicer field
EeffParamLim1_non = 1.e2;  % Lower limit
EeffParamLim2_non= 1.e4;   % Upper limit

% Normal distribution with  (1) mu=Eeff_mean (2) ?=Eeef_sigma
Eeff_normal=false;         % then Eeff from uniform distrib

Leff_plaw=false;           % If true Leff : from power law
                           % if false Leff: as a linear function of Eeff
                          
Eeff_mean_coeff=1.e4;     % Multiplied by ED to give Eeff_mean
Eeef_sigma_coeff=(1/800)*Eeff_mean_coeff; % Multiplied by ED to give Eeff_mean

check_EeffArr1=false;      % to check Eeff, Leff,Beff

% Uniform Distribution: Eeff = ED*Efactor_non
Efactor_non = 1.e5;        % When constant, Eeff = (Dreicer field) x Efactor []

if Eeff_normal  
 lEplaw_non=false;         % We can't have both 
end
%% OUR MODEL DOES NOT TAKE INTO ACCOUNT COLLISIONS
% Collisions
lbackIon = false;       % Background species: ions (T) or electrons (F)
TCollFactor = 1;        % Backround species temperature (in units of initial temperature)

% Timing
lTime = true;           % Measure execution time (for-loop)

% Debugging
lstats = false;
lCstats = false;
lCSstats = false;
if check_EeffArr1
%BeffArr1=zeros(1,2*10^5);
EeffArr1=zeros(1,2*10^5);
LeffArr1=zeros(1,2*10^5);  %dont forget j1=1; under ip=1:nP
end
%}
%--------------------------------------------------------------------------
%N = N + 1 - mod(N,2);   % ensure that N is odd
gammaLim = 1.00001;     % Relativistic/classical limit

%% ===== END OF USER DEFINED PARAMETERS
%% ========================================================================
%% Constants
me = 9.10938e-28;       % rest mass of the electron [gr]
mp = 1.67262e-24;       % rest mass of the proton [gr]
qe = 4.803204e-10;      % elementary charge [statC]
c  = 2.99792458e10;     % speed of light [cm/sec]
kB = 1.38064e-16;       % Boltzmann's constant [erg/K]
erg2eV = 6.2415e11;     % erg to eV conversion factor [eV/erg]
% In SI units
q_SI = 1.6022e-19;      % elementary charge [C]
kB_SI = 1.3807e-23;     % Boltzmann constant [J/K]
eps0 = 8.8542e-12;      % Dielectric constant of free space [F/m]
%% ========================================================================
%% ===== Derived parameters and switches
VA    = B/sqrt(4*pi*n0*mp);                 % Alfven speed [cm/sec]
vth_e = sqrt(2 * kT0 ./erg2eV ./ me);       % electron thermal velocity [cm/sec]
vth_p = sqrt(2 * kT0 ./erg2eV ./ mp);       % proton thermal velocity [cm/sec]
vth_i = sqrt(2 * kT0 ./erg2eV ./ (Zion*mp));% ion thermal velocity [cm/sec]
c2    = c.^2;

if lIon
    m = Zion * mp;          % particle mass [gr]
    q = Zion * qe;          % particle charge [statC]
    vth = vth_i;            % particle thermal velocity [cm/sec]
else
    m = me;                 % particle mass [gr]
    q = -qe;                % particle charge [statC]
    vth = vth_e;            % particle thermal velocity [cm/sec]
end

Wrest = m*c2;              % particle rest energy [erg]
gamma = 1 ./ sqrt(1-(vth./c).^2);
if gamma > gammaLim             % initial energy distribution [erg]
    Wth = (gamma-1) .* Wrest;
else
    Wth = m*vth.^2 / 2;
end

V = Vfactor * VA;                  % speed of the cloud [cm/sec]
GammaV = 1 ./ sqrt(1-V.^2 ./c2);   % Lorentz factor for the cloud

% Scattering switch
lscattSwitch = true;
if imodel == 0
    lscattSwitch = false;
end
% Collisions switch
lcollSwitch = false;
if iCollFlag ~= 0
    lcollSwitch = true;
end

%% ------------------------------------------------------------------------
%% ----- Grid   # NOT USED 
%   Used for directioning
twondims = 2*ndims;
dimsArr = [-ndims:-1 1:ndims];
%% ------------------------------------------------------------------------
%% ----- Particles

% === Initial energy
if lW0const     % constant (monoenergetic)
    W0_const = kT0 ./ erg2eV;
    gamma    = 1 + W0_const ./ Wrest;
    if gamma > gammaLim                 % initial velocity [cm/sec]
        v0_const = c * sqrt(1 - 1 / gamma^2);
    else
        v0_const = sqrt(2*W0_const / m);
    end
    W0 = W0_const * ones(1,nP);         % initial energy array [erg]
    v0 = v0_const * ones(1,nP);         % initial velocity array [cm/sec] %*** Needed?
else            % maxwellian
    mu = 0;
    sig2 = (kT0 ./ erg2eV) / m;         % variance; first convert eV to erg
    v0 = abs(mu + sqrt(sig2) * randn(1,nP));  % initial velocity distribution [cm/sec]

    % Apply a threshold to the initial speed
    [nn, ~] = find(abs(v0(:)) < v0thresFactor * vth);
    for ii = 1 : length(nn)
        while abs(v0(nn(ii))) < v0thresFactor * vth
            v0(nn(ii)) = abs( mu + sqrt(sig2) * randn );
        end
    end
   
    W0(nP) = 0;
    for ip = 1:nP
        gamma = 1 ./ sqrt(1-(v0(ip)./c).^2);
        if gamma > gammaLim             % initial energy distribution [erg]
            W0(ip) = (gamma-1) .* Wrest;
        else
            W0(ip) = m*v0(ip).^2 / 2;
        end
    end
end
%% ------------------------------------------------------------------------
%% ----- Collisions -> NOT USED
% Coulomb logarithm; used for collisions and Dreicer field
ne = n0;
ni = n0;
n1 = n0;
if lbackIon
    m2 = Zion * mp;
    n2 = ni;
else
    m2 = me;
    n2 = ne;
end
E1_SI = kT0 ./ erg2eV * 1.e-7;  % Foreground species energy [J]
n1_SI = n1 * 1.e6;              % Foreground species density [m^-3]
n2_SI = n2 * 1.e6;              % Background species density [m^-3]
m_SI  = m * 1.e-3;              % Foreground species mass [kg]
Temp1 = E1_SI ./ kB_SI;         % Foreground species temperature [K]
Temp2 = Temp1 .* TCollFactor;   % Background species temperature [K]
E2_SI = Temp2 .* kB_SI;         % Background species energy
% Coulomb logarithm for collisions of electrons over ions.
%   It is similar for either e/e or i/i of the same energy.
%*** Zion???
lnCoul = log(12*pi * (eps0 .* kB_SI * Temp2)^1.5 ./ (q_SI^3 .* sqrt(n2_SI)));
% This is the Coulomb logarithm for collisions of electrons over ions.
% However, it is similar for either e/e or i/i of the same energy

% Other collision parameters
if lcollSwitch
    % SI units
    Gcoll12_SI = n2_SI .* q_SI.^4 .* lnCoul ./ (4*pi .* eps0.^2 .* m_SI.^2);
    vth_SI = sqrt(2 * kB_SI * Temp1 / m_SI);
    % Even though derived from SI formulas, the following have units sec^-1
    % and sec, thus suitable for CGS.
    nuColl = Gcoll12_SI ./ vth_SI^3;    % Collision frequency
    tColl = 1 / nuColl;                 % Collision time

    % Return to CGS
    mfpColl = tColl * vth;              % Collision mean free path [cm]
    kTm = kB * Temp2 / m;
    kstepColl(nP) = 0;                  % Number of collision intervals per particle
    nuColl_arr(nP) = 0;                 % Sum of nu per particle (to estimate a mean)

    if lCstats
        iter1Arr = [];
        iter2Arr = [];
        tau1Arr = [];
        tau2Arr = [];
        tauArr = [];
        root1Arr = [];
        root2Arr = [];
        rootsArr = [];
        nuCollArr = [];
        lFreeArr = [];
        muArr = [];
        vtauArr = [];
        dvtauArr = [];
    end
end
clear ne ni n1 n2 m2 n1_SI n2_SI m_SI E1_SI E2_SI Temp1 vth_SI
%% ------------------------------------------------------------------------
%% ----- Reconnecting Current sheets
% Dreicer field
ED = 2.33e-8 * (1.e-9 * n0) * (1.e7 / T0cs) * (lnCoul / 23.2);   % NS for dreicer field

if lBdepend
    alpha1=c*(LeffLim2-LeffLim1)/(VA*(BLim2-BLim1));
    beta1=LeffLim1-BLim1*((LeffLim2-LeffLim1)/(BLim2-BLim1)); 
   
    if lCSstats
      BeffArr = [];
       EeffArr = [];
       LeffArr = [];
    end
end
if lstats
    nScat1(nP) = 0;     % ***debug
    ngamma(nP) = 0;     % relativistic scattering
    maxgamma(nP) = 0;   % ... statistics
    gammaArr = [];
    ndtScat(nP) = 0;    % timed scattering
    maxdtScat(nP) = 0;  % ... statistics
    dtCSArr = [];
end
%% -----NON Reconnecting Current sheets
    
     if lEplaw_non      % Power-Law distribution
            EeffLim1_non = EeffParamLim1_non * ED;
            EeffLim2_non = EeffParamLim2_non * ED;    
            alpha_non=(LeffLim2_non-LeffLim1_non)/(EeffLim2_non-EeffLim1_non);
            beta_non=LeffLim1_non-EeffLim1_non*((LeffLim2_non-LeffLim1_non)/(EeffLim2_non-EeffLim1_non));
        if lCSstats, EeffArr = []; end
     elseif Eeff_normal  %Normal Distribution
         Eeff_mean=Eeff_mean_coeff*ED;
         Eeef_sigma=Eeef_sigma_coeff*ED;
         
         alpha_non1=(LeffLim2_non-LeffLim1_non)/((Eeff_mean+3*Eeef_sigma)-(Eeff_mean-3*Eeef_sigma));
%         beta_non1=LeffLim1_non-EeffLim1_non*((LeffLim2_non-LeffLim1_non)/(EeffLim2_non-EeffLim1_non));
         beta_non1=LeffLim1_non-(Eeff_mean-3*Eeef_sigma)*((LeffLim2_non-LeffLim1_non)/((Eeff_mean+3*Eeef_sigma)-(Eeff_mean-3*Eeef_sigma)));
     
     else                % Uniform distribution
         
        Eeff = Efactor_non * ED;
    end
    if lLeffrand_non
        if lCSstats, LeffArr = []; end
    else
        Leff = LeffLim1_non;
    end

if lstats
    nScat1(nP) = 0;     % ***debug
    ngamma(nP) = 0;     % relativistic scattering
    maxgamma(nP) = 0;   % ... statistics
    gammaArr = [];
    ndtScat(nP) = 0;    % timed scattering
    maxdtScat(nP) = 0;  % ... statistics
    dtCSArr = [];
end
%% ------------------------------------------------------------------------
%% ----- Statistics and Monitoring
% ==== Initialize monitoring
lmonitSwitch = true;
if ntM == 0 || abs(tfin-tinit) < 1.e-10
    lmonitSwitch = false;
    ltfinStop = false;
end
if lmonitSwitch
    dtM = (tfin-tinit) / ntM;   % monitoring time step
    tM = dtM:dtM:tfin;          % monitoring time array
    %tM=linspace(t_init_monit,tfin,ntM); %                [NS]
    tM=logspace(log10(t_init_monit),log10(tfin),ntM);
    %kpM(ntM) = 0;               % particles present in each monitoring
    rmsd(nP,ntM) = 0;           % Mean squared displacement per particle using euclidean distance
    %rmsd_(nP,ntM) = 0;          % Mean squared displacement per particle using euclidean distance
    %rmsd1(ntM) = 0;             % Mean squared displacement per particle using euclidean distance
    Wkinet(nP,ntM) = 0;         % Kinetic energy of the particles
    %Vkinet(nP,ntM) = 0;          % NS
    %ktkicks(nP,ntM) = 0;        % Number of kicks, i.e. visits to an ACTIVE...
                                %   grid point (acceleration event)
     %rdist3(nP,ntM)=0;                         
   % part2rms(nP,ntM,3) = 0; %HI
   % is_in = zeros(nP,ntM);  %HI: this is a marker, set to 1
                            %    if a particle is in the system
                            %    at a given monitoring time,
                            %    else it is 0
                            
%       a11=zeros(nP,20*N1);                     
end

%disp(ntM)
%disp(size(rmsd))

% ==== Escape
%t_esc_new(nP)=0;            % Escape time for particles leavig the box
t_esc(nP) = 0;              % Escape time per particle
W_esc(nP) = 0;              % Kinetic energy of the particle on escape
%v_esc(nP) = 0;              % Velocity of the particle on escape

% ==== Statistics
kkicks(nP) = 0;             % number of kicks, i.e. visits to an ACTIVE...
                            %   grid point (acceleration event)
%forwkicks(nP) = 0;          % Counter for following collisions per particle
%hdonkicks(nP) = 0;          % Counter for head on collisions per particle

%kFreePath(nP) = 0;          % mean free path per particle
mDelW_W(nP) = 0;            % < (\Delta W) / W > per particle

sumdW_p(nP) = 0;            % Sum of dW per particle
%countdW_p(nP) = 0;

%nSignore(nP) = 0;           % total number of ignored scatterings
                            %                       and reinserted particles
%iSignore = uint32([]);      % indices number of reinserted particles;
                            %           0 -> (2^32-1) = 4.294.967.295
%tSignore = [];              % times of ignored scatterings

if lcollSwitch
    nCignore(nP) = 0;       % total number of reinserted particles
    iCignore = uint32([]);  % indices number of reinserted particles;
                            %           0 -> (2^32-1) = 4.294.967.295
    tCignore = [];          % exit time of reinserted particles
    rootInaccur = [];
end

%% Complete timeseries
%if lTScomplete
   % W_Cell = {};
   % t_Cell = {};
    dW_Cell = {};
   % dW_Cell_overW= {};
   % dtScat_Cell = {};
   % kicks_Cell = {};
   % if imodel > 10, model_Cell = {}; end
%end

%% ------------------------------------------------------------------------
%% ----- Export setup parameters

%{


%NAct = nnz(isActive);       % number of active points

%mfp0 = l./ActiveRatio;      % theoretical mean free path [cm]
%alpha0 = 4/3*V^2/(c*mfp0);  %
%tacc0 = 1/alpha0;           % theoretical acceleration time [sec]

model.imodel = imodel;
model.ndims = ndims;
model.iBCs = iBCs;
model.iCollFlag = iCollFlag;
model.const.gammaLim = gammaLim;
model.const.erg2eV = erg2eV;

scatterers.B = B;
scatterers.N = N;
scatterers.L = L;
%scatterers.R = ActiveRatio;
scatterers.n0 = n0;
scatterers.Vfactor = Vfactor;
scatterers.VA = VA;
scatterers.l = l;
%scatterers.NAct = NAct;
%scatterers.mfp0 = mfp0;
%scatterers.tacc0 = tacc0;

particles.nP = nP;
particles.kT0 = kT0;
particles.vth = vth;
particles.Wth = Wth;
particles.Wrest = Wrest;
particles.lIon = lIon;
particles.Zion = Zion;
particles.l0centre = l0centre;
particles.lW0const = lW0const;
particles.v0thresFactor = v0thresFactor;

monitor.lmonitSwitch = lmonitSwitch;
if lmonitSwitch
    monitor.lmonitSwitch = lmonitSwitch;
    monitor.tinit = tinit;
    monitor.tfin = tfin;
    monitor.ntM = ntM;
    monitor.WmonitFactor = WmonitFactor;
end

collisions.iCollFlag = iCollFlag;
if lcollSwitch
    collisions.lbackIon = lbackIon;
    collisions.TCollFactor = TCollFactor;
    collisions.tColl = tColl;
    collisions.mfpColl = mfpColl;
    collisions.lnCoul = lnCoul;
end

if imodel == 2 || imodel == 12
    EF.idtCS = idtCS;
    EF.T0cs = T0cs;
    EF.ED = ED;
    EF.lBdepend = lBdepend;
    EF.lEplaw = lEplaw;
    EF.lLeffrand = lLeffrand;
    EF.LeffLim1 = LeffLim1;
    EF.LeffLim2 = LeffLim2;
    EF.BPLa = BPLa;
    EF.lBPLparam = lEPLparam;
    EF.BLim1 = BLim1;
    EF.BLim2 = BLim2;
    EF.EPLa = EPLa;
    EF.lEPLparam = lEPLparam;
    EF.EeffParamLim1 = EeffParamLim1;
    EF.EeffParamLim2 = EeffParamLim2;
    if ~lBdepend && lEplaw
        EF.EeffLim1 = EeffLim1;
        EF.EeffLim2 = EeffLim2;
    end
    EF.Efactor = Efactor;
    model.UCS = EF;
elseif imodel == 3 || imodel == 13
    CD.ldWTerm_cosPhi = ldWTerm_cosPhi;
    CD.ldWTerm_dBB = ldWTerm_dBB;
    CD.idWTerm_W = idWTerm_W;
    CD.idWTerm_const = idWTerm_const;
    CD.BPLa = BPLa;
    CD.lBPLparam = lEPLparam;
    CD.BLim1 = BLim1;
    CD.BLim2 = BLim2;
    model.UCS = CD;
end

%save('setup006.mat',...
 %  'model', 'particles', 'scatterers', 'collisions', 'monitor'...
  %);

%}
% Initial data
W0_eV = W0 .* erg2eV;
%save('init006.mat', 'W0_eV', 'v0');
%keyboard
clear  alpha0 ...
        model particles scatterers collisions monitor %mfp0  tacc0 NAct
%--------------------------------------------------------------------------
%% ========================================================================
%% ===== MAIN
%% Particle loop
disp('Starting particle loop')
if lTime, tic; toc0 = toc; end

    %{
N_check=0; % for the case where number of kicks is equal to N
    N_check2=0;           % we re-run the simulation for that particle with a new
                 % N=2N
      %}           
    j1=1; %Only if you want to check the values for Beff, Leff, Eeff (model 2)            
for ip = 1:nP
      
         N=N1; % the initial value of N       %NS
         
        Dw_at_t0_is_0= Dw_at_t0_is_0_or_Dw; % as said above not enrg at t=0
 
    if mod(ip,1000) == 0  % Just to check progress
        if lTime
            toc1 = toc;
           disp(['Particle: ' num2str(ip) '; elapsed time [sec]: ' num2str(toc1-toc0) ', total: ' num2str(toc1)])
            toc0 = toc1;
        else
            %disp(['Particle: ' num2str(ip)])
        end
    end

    
%% Prelocate Particle coordinates -> Important for MONTE CARLO
x = zeros(1,N);                            % Prelocate coordinates
y = zeros(1,N);                            % Prelocate coordinates              NS
z = zeros(1,N);                            % Prelocate coordinates
costheta = -1 + (1-(-1)).* rand(1,N-1);    %Generate random numbers between -1 < cos(theta) < 1   NS
phi      = 360*rand(1,N-1);                %Generate random angle phi 0< phi < 2*Pi     NS
theta    = (180/pi).*acos(costheta);  
sintheta = sind(theta);
cosphi   = cosd(phi);                  
sinphi   = sind(phi);

 if l0centre  
   x(1,1)=0;
   y(1,1)=0;     % Select  initial partice position at (0,0,0) NS
   z(1,1)=0;
 else
  
   x(1,1) = randi([-L2,L2],1,1);
   y(1,1) = randi([-L2,L2],1,1);   % Select  initial partice position randomly 

    if square_to_rectangle
       z(1,1) = randi([-L2z,L2z],1,1);
    else
       z(1,1) = randi([-L2,L2],1,1);
    end
 end
%% Just to check

%tau11=zeros(1,N);   % see line 1091, used to check the time between the last energ and the escape
 
tau_check=zeros(1,N);
%speed_check=zeros(nP,N,3);


%% Generating power-law distributed random numbers with the method of transformation
 dim=[1,N-1];        % dimensions of dr matrix
   
 % The steps of each particle are dictated by the distribution
 % P(dr) ~ dr^(P_law_index). The following script produces those steps
 if upper_step_limit 

   dr=((Dr_min^(1-P_law_index)) + (Dr_max^(1-P_law_index)-Dr_min^(1-P_law_index)).*rand(dim)).^(1/(1-P_law_index));
 else
   dr=Dr_min.*(1-rand(dim)).^(1/(1-P_law_index));
 end

%% Compute coordinates of particles for each iteration
 for j=2:N
   x(1,j)=x(1,j-1)+dr(1,j-1)*sintheta(1,j-1)*cosphi(1,j-1);
   y(1,j)=y(1,j-1)+dr(1,j-1)*sintheta(1,j-1)*sinphi(1,j-1);         % Compute coordinates of particles for each iteration NS
   z(1,j)=z(1,j-1)+dr(1,j-1)*costheta(1,j-1);  
 end 
    %% Initialize particle position
    %# Initialize particle energy, velocity, direction and time
    W = W0(ip);
    W_prev = W; % Same cmd
    v = v0(ip);                                                             
    % direction % Next two lines ->Not needed (remnants of Latice model)
    pdir = zeros(1,ndims); icoin = dimsArr(randi(twondims,1));
    pdir(abs(icoin)) = sign(icoin);

   % Important for monitoring the particles at predifined monitoring times
    t = tinit;
    t_prev = t;
    t_prev_M = t;   %HI: note the difference between T_prev and t_prev_M

    kFree = 1;      % number of free steps since last acceleration
    DelW_W = 0;     % \sum (\Delta W) / W

    % Initialize monitoring
    if lmonitSwitch
        xy(1) = x(1,1);                %NS
        xy(2) = y(1,1);
        if (ndims == 3), xy(3) =z(1,1); end
    end

    %***
    sumdW_p(ip) = 0;
    %countdW_p(ip) = 0;

    % Complete timeseries
    %if lTScomplete
        %W_Arr = W;
        %t_Arr = tinit;
        dW_Arr = NaN;
        %dW_Arr_overW = NaN;
        %dtScat_Arr = NaN;
       % kicks_Arr = 0;
       % if imodel > 10, model_Arr = NaN; end
   % end
    
    %% particle travels in the domain
   for j=(2-energ_t0):N   % if energ_t0=1 (energ at t=0) NS
       
      if j>energ_t0       %so at t=0 I have an energization event
        if t >= t_final    %    tfin
            if iBCs == 1 || iBCs == 2, break; end              

            if iBCs == 0 && ltfinStop,  break; end
        end
     %% Particle reached the boundaries
   if ((x(1,j-1)<=dmin || x(1,j-1)>=dmax || y(1,j-1)<=dmin || y(1,j-1)>=dmax) ||...
                (ndims==3 && (z(1,j-1)<=dmin_z || z(1,j-1)>=dmax_z)))
            if iBCs == 0        % Stop
                %disp(num2str(j))
                break
            elseif iBCs == 1    % Enter from the oposite  ( We never used periodic boundary conditions
                if (i==1 && pdir(1)==-1), i = N; end      % So i never had to change from latice to MOnte carlo
                if (i==N && pdir(1)==1) , i = 1; end      %  for periodic boundary conditions)
                if (j==1 && pdir(2)==-1), j = N; end
                    if (j==N && pdir(2)==1) , j = 1; end
                if (ndims == 3)
                    if (k==1 && pdir(3)==-1), k = N; end
                    if (k==N && pdir(3)==1) , k = 1; end
                end
            elseif iBCs == 2    % New position, new velocity
                i = randi([2 N-1],1); % New position not at the boundary
                j = randi([2 N-1],1);
                if (ndims == 3), k = randi([2 N-1],1); end
                % New energy from maxwellian
                v = abs(sqrt(sig2) * randn);        % new velocity [cm/sec]
                while abs(v) < v0thresFactor * vth  % apply a threshold
                    v = abs(sqrt(sig2) * randn);
                end
               
                gamma = 1 ./ sqrt(1-(v./c).^2);
                if gamma > gammaLim
                    W = (gamma-1) .* Wrest;
                else
                    W = m*v.^2 / 2;
                end
            end
    end     % EO boundary check
       
    %% update position
        
        xy_prev = xy ;    %HI: needed for linear interpolation
        if lmonitSwitch   % Particle position; used only for r-distance
            if iBCs == 0
                xy(1) = x(1,j);
                xy(2) = y(1,j);   %NS
                if (ndims == 3), xy(3) = z(1,j); end
            else
                xy(1) = xy(1) + l * pdir(1);              
                xy(2) = xy(2) + l * pdir(2);
                if (ndims == 3), xy(3) = xy(3) + l * pdir(3); end
            end
        end

        W_prev = W;
        t_prev = t;

        %% Active-Inactive check
      %  if (isActive(i,j,k)) % && abs(v) >= ScatLim)
      %  ?? need to check always active point    NS
            %% collisions
            if lcollSwitch
                lFree = l * kFree;
                kstepColl(ip) = kstepColl(ip) + 1;

                N1 = randn;
                N2 = abs(randn) * sign(v);
                %clear rootCollFUN
                nuColl_v = nuColl .* vth.^3 ./ abs(v).^3;
                nuColl_arr(ip) = nuColl_arr(ip) + nuColl_v;
                rootCollFUN = ...
                    @(tau) abs(StauFUN(tau, v, nuColl_v, N1, N2, kTm)) - lFree;
                tau1 = lFree / abs(v);
                tau2 = sqrt(lFree / abs(sqrt(kTm)*N1));
                if tau1 > tau2;
                    tau12 = tau1;
                    tau1 = tau2;
                    tau2 = tau12;
                end
                iter1 = 0;
                iter2 = 0;
                root1 = rootCollFUN(tau1); % Fewer function calls in 2nd while
                while root1 > 0
                    tau1 = 0.1 * tau1;
                    iter1 = iter1 + 1;
                    root1 = rootCollFUN(tau1);
                end
                root2 = rootCollFUN(tau2);
                while root2 < 0
                    %rootCollFUN(tau2)
                    tau2 = 2 * tau2;
                    iter2 = iter2 + 1;
                    root2 = rootCollFUN(tau2);
                end
                tau = fzero(rootCollFUN,[tau1 tau2]);
                mu = exp(-nuColl_v * tau);
                vtau = v * mu + sqrt(kTm * (1-mu^2)) * N1;

                if lCstats
                    iter1Arr = [iter1Arr; iter1];
                    iter2Arr = [iter2Arr; iter2];
                    tau1Arr = [tau1Arr; tau1];
                    tau2Arr = [tau2Arr; tau2];
                    root1Arr = [root1Arr; root1];
                    root2Arr = [root2Arr; root2];
                    tauArr = [tauArr; tau];
                    rootsArr = [rootsArr; rootCollFUN(tau)];
                    lFreeArr = [lFreeArr; lFree];
                    vtauArr = [vtauArr; vtau];
                    dvtauArr = [dvtauArr; vtau-v];
                    nuCollArr = [nuCollArr; nuColl_v];
                    muArr = [muArr; mu];
                end
                %v = abs(vtau);
                v = vtau;
                gamma = 1 ./ sqrt(1-(v./c).^2);
                if gamma > gammaLim
                    W = (gamma-1) .* Wrest;
                else
                    W = m*v.^2 / 2;
                end
               
                if gamma < 1    % Particle velocity is too low
                    W = W_prev;
                    nCignore(ip) = nCignore(ip) + 1;
                    iCignore = [iCignore; ip];
                    tCignore = [tCignore; t];
                end
                %t = t + tau;

                if abs(rootCollFUN(tau)) > 1.e-6
                    rootInaccur = [rootInaccur; ip rootCollFUN(tau)];
                end
            else
               tau = abs(dr(1,j-1))/abs(v);  %NS
               %v_check(ip,j-1)=v;
                tau_check(1,j-1)=tau;
            end % EO collisions
            t = t + tau;
            %t_mat(ip,j)=t;
            
     end         
            %% scattering event  % energization of particle
            if lscattSwitch
                % Scattering event
                kkicks(ip) = kkicks(ip) + 1;
              
                % Scatterer direction
                sdir = zeros(1,ndims); icoin = dimsArr(randi(twondims,1));
                sdir(abs(icoin)) = sign(icoin);
                
               
                isgn = pdir(1)*sdir(1) + pdir(2)*sdir(2); % Faster than dot(pdir,sdir)
                if (ndims == 3), isgn = isgn + pdir(3)*sdir(3); end
                
                      % end                        
%% update Possibility to encounter ASs  
        if magnetic_field_decay  
          PdW=PdW_init*log10(t/t_decay);     % decaying magnetic field       NS  
        end
%%
                Wi = W;             % Kinetic
                WiTot = Wi + Wrest; % Total

                if Dw_at_t0_is_0>0
                
                % Energy and time change
                if imodel == 12
                    idWChoice =     (rand>PdW) + 1; % (0,1) -> (1,2)
                elseif imodel == 13
                    idWChoice = 2 * (rand>PdW) + 1; % (0,1) -> (1,3)
                else
                    idWChoice = imodel;
                end
                   
                switch idWChoice
                    case 1
                        dW = 2*GammaV * WiTot .* (V^2./c2 - isgn * (V*v)./c2);
                        
                        %dW=0;
                        dtScat = 0;
                    case 2
                                                        %NS for
                    if  rand >P_reconnect   % Non-recommecting UCS
                       
                        if lEplaw_non   %power law distribution 
                             Eeff =((EeffLim1_non^(1-EPLa_non)) + (EeffLim2_non^(1-EPLa_non)-EeffLim1_non^(1-EPLa_non)).*rand).^(1/(1-EPLa_non)); %NS
                             Leff = Eeff*(alpha_non)+beta_non;                                                                  %NS
                             sign1 =(rand(1,1) > 0.5)*2 - 1;
                             
                      elseif Eeff_normal   % Normal Distribution
                                   
                                   Eeff=Eeff_mean-3*Eeef_sigma-1;
                                   while (Eeff<Eeff_mean-3*Eeef_sigma) || (Eeff > Eeff_mean+3*Eeef_sigma)
                                     Eeff = normrnd(Eeff_mean,Eeef_sigma,1); 
                                   end
                                   %Eeff = Eeff_mean+TruncatedGaussian(Eeef_sigma,[Eeff_mean-3*Eeef_sigma, Eeff_mean+3*Eeef_sigma]-Eeff_mean,1);
                                  if Leff_plaw
                                     Leff =((LeffLim1_non^(1-EPLa_non)) + (LeffLim2_non^(1-EPLa_non)-LeffLim1_non^(1-EPLa_non)).*rand).^(1/(1-EPLa_non));  % Power law   
                                  else
                                     Leff = Eeff*(alpha_non1)+beta_non1;  
                                  end
                                  sign1 =(rand(1,1) > 0.5)*2 - 1;
                       
                        else       % Uniform distribution
                                   
                                  Eeff = Efactor_non * ED; 
                                  Leff =((LeffLim1_non^(1-EPLa_non)) + (LeffLim2_non^(1-EPLa_non)-LeffLim1_non^(1-EPLa_non)).*rand).^(1/(1-EPLa_non));  % Power law   
                                  sign1 =(rand(1,1) > 0.5)*2 - 1; 
                         end
                               
                         
                    else        % Reconnecting UCS
          
                            Beff=((BLim1^(1-BPLa)) + (BLim2^(1-BPLa)-BLim1^(1-BPLa)).*rand).^(1/(1-BPLa)); % NS
                            Eeff = Beff .* VA ./ c;
                            Leff = Eeff *alpha1 + beta1;
                            sign1=1;
                          
                           % BeffArr1(j1,ip)=Beff;
                            if lCSstats
                                BeffArr = [BeffArr; Beff];
                            end
                       
                    end
                     
                            if lCSstats
                                 EeffArr = [EeffArr; Eeff];
                                LeffArr = [LeffArr; Leff];
                            end
                            
                        dW = sign1.*abs(q) .* Eeff .* Leff;
                        dtScat = idtCS * Leff ./ sqrt(2*Wi/m);
                        
                        if lstats
                            nScat1(ip) = nScat1(ip) + 1;
                            ndtScat(ip) = ndtScat(ip) + 1;
                            maxdtScat(ip) = max(maxdtScat(ip), dtScat);
                            dtCSArr = [dtCSArr; dtScat];
                            % Check if relativistic
                            gamma = 1 + Wi/Wrest;
                            if gamma > gammaLim
                                ngamma(ip) = ngamma(ip) + 1;
                                maxgamma(ip) = max(maxgamma(ip), gamma);
                                gammaArr = [gammaArr; gamma];
                            end
                        end
                    case 3
                        randTerm = 1;
                        if ldWTerm_cosPhi
                            randTerm = randTerm .* cos(pi * rand - pi/2);
                        end
                        if ldWTerm_dBB
                            dB = (BLim1^(-BPLa) - BPLa/BPLnorm *rand).^(-1/BPLa);
                            randTerm = randTerm .* dB ./ B;
                        end
                       
                        dW = 2*GammaV * (idWTerm_W .* WiTot + (idWTerm_W-1) * m.*c2) .*...
                                           (V.^2./c2 .* idWTerm_const  + ...
                                            V.*v./c2 .* randTerm);
                        dtScat = 0;%idtCS * Leff ./ vi;
                end
               
                          if check_EeffArr1
                         EeffArr1(j1)=Eeff;
                          LeffArr1(j1)=Leff;
                          j1=j1+1;
                          end
                
                
                % Complete timeseries
               % if lTScomplete
                   % W_Arr = [W_Arr Wi];
                    dW_Arr = [dW_Arr dW];
                    %dW_Arr_overW = [dW_Arr_overW a11];
                    
                    %t_Arr = [t_Arr t];
                   % dtScat_Arr = [dtScat_Arr dtScat];
                   % kicks_Arr = [kicks_Arr kkicks(ip)];
                   % if imodel > 10, model_Arr = [model_Arr idWChoice]; end
               % end
               
                   Wf = Wi + dW;   % Kinetic energy
                   %t = t + dtScat; 
               %    a11(ip,j)=abs(dW/WiTot);
                   
                    
                   if (Wf<0) && (second_way)
                                                         % NS: For non-rec UCS
                  Wf=abs(Wf);                            % sometimes Dw'<0
                   end
                     
               
                t = t + dtScat; %
                %t_mat(ip,j-1)=t;
               
                sumdW_p(ip) = sumdW_p(ip) + dW;
                %countdW_p(ip) = countdW_p(ip) + 1;
                % Calculate particle velocity
                gamma = 1 + Wf/Wrest;
               
                if gamma < 1    % Particle velocity is too low
                    % Instead of reinserting, ignore last collision and continue
                    Wf = Wi;    % Kinetic
                   % nSignore(ip) = nSignore(ip) + 1;
                    %iSignore = [iSignore; ip];
                    %tSignore = [tSignore; t];
                   
                    % Complete timeseries
                   % if lTScomplete
                      %  W_Arr = W_Arr(1:end-1);
                        dW_Arr = dW_Arr(1:end-1);
                       % dW_Arr_overW=dW_Arr_overW(1:end-1);
                        %t_Arr = t_Arr(1:end-1);
                        %dtScat_Arr = dtScat_Arr(1:end-1);
                        %kicks_Arr = kicks_Arr(1:end-1);
                        %if imodel > 10, model_Arr = model_Arr(1:end-1); end
                    %end
                end
               
                % Counters for statistics
                if imodel == 1
                    if isgn ~= 0
                        if isgn > 0
                           % forwkicks(ip) = forwkicks(ip) + 1 ;
                        else
                           % hdonkicks(ip) = hdonkicks(ip) + 1 ;
                        end
                    end
                end
                %kFreePath(ip) = kFreePath(ip) + kFree;
                kFree = 1;
                DelW_W = DelW_W + abs(dW)/ WiTot;

                % New energy, velocity and direction
                W = Wf;     % Kinetic energy
                if gamma > gammaLim
                    v = c * sqrt(1 - 1 / gamma^2);
                else
                    v = sqrt(2*W/m);
                end
                
               %{ 
                if j==1
                    if idWChoice==2
                        dtScat =  Leff ./ sqrt(2*Wi/m);
                    else
                        L =((LeffLim1_non^(1-EPLa_non)) + (LeffLim2_non^(1-EPLa_non)-LeffLim1_non^(1-EPLa_non)).*rand).^(1/(1-EPLa_non));
                    dtScat= abs(L)/abs(v-v0(ip));  %NS
                    end
                    t = t + dtScat;
                end
               %}
                pdir = zeros(1,ndims); icoin = dimsArr(randi(twondims,1));
                pdir(abs(icoin)) = sign(icoin);
                
                 end
                Dw_at_t0_is_0=1;
            end % EO scattering

        %% Monitoring
               
        if lmonitSwitch && (W >= WmonitFactor * Wth) % Apply a threshold to the energy
                                                     % of the particles being monitored
            % Find correct monitor time
            [~,itM] = min(abs(t-tM));
            [~,itM_prev] = min(abs(t_prev_M-tM)); %HI t_prev_M, not t_prev
            if (tM(itM) > t)
                itM = itM - 1;
            end
            if (tM(itM_prev) < t_prev_M) %HI t_prev_M
                itM_prev = itM_prev + 1;
            end

            % Update monitor time, if necessary
            %disp(size(rmsd))
            if (itM>=itM_prev)
                %HI linear interpolation of position
                speed1 = (xy(1) - xy_prev(1))/(t-t_prev) ; %HI
                speed2 = (xy(2) - xy_prev(2))/(t-t_prev) ; %HI
                speed3 = (xy(3) - xy_prev(3))/(t-t_prev) ; %HI
                
              
                for jtM = itM_prev : itM
                    xy1 = xy_prev(1) + speed1*(tM(jtM) - t_prev); %HI
                    xy2 = xy_prev(2) + speed2*(tM(jtM) - t_prev); %HI
                    xy3 = xy_prev(3) + speed3*(tM(jtM) - t_prev); %HI
                    %speeda(ip,jtM)=speed1;
                   % speedb(ip,jtM)=speed2;
                   % speedc(ip,jtM)=speed3;
                    if (tM(jtM)-t_prev < 0)
                        disp(tM(jtM)-t_prev)
                    end
                    if (abs(xy1 - xy_prev(1)) < 1.0e-10 ...
                      & abs(xy2 - xy_prev(2)) < 1.0e-10...
                      & abs(xy3 - xy_prev(3)) < 1.0e-10)
                        disp(xy1)
                        disp(xy1-xy_prev(1))
                    end

                    t_prev_M = t ; %HI: different from t_prev !
                    
                    if ((xy1<=dmin || xy1>=dmax || xy2<=dmin || xy2>=dmax) ||...
                         (ndims==3 && xy3<=dmin_z || xy3>=dmax_z))
                     
                     rmsd(ip,jtM) = 0;                              % Space msd
                    %rmsd_(ip,jtM) = 0 ;          % Space msd
                    %rmsd1(jtM) = 0; % Space msd
       
                    Wkinet(ip,jtM) = 0;         % Evolution of kinetic energy
                    %is_in(ip,jtM) = 0;                                         %NS: the particle is outside
                        
                     
                    else
                          
                    %kpM(jtM) = kpM(jtM) + 1;    % Number of particles in this time
                   % ktkicks(ip,jtM) = kkicks(ip);
                    rdist2 = (xy1 - x(1,1)).^2 + (xy2-y(1,1)).^2;   % Eucledian           %NS
                     
                      if (ndims == 3), rdist2 = rdist2 + (xy3 - z(1,1)).^2; end 
                    rmsd(ip,jtM) = rdist2;                              % Space msd
                    %rmsd_(ip,jtM) = rmsd_(ip,jtM) + rdist2 ;            % Space msd
                  %  rmsd1(jtM) = rmsd1(jtM) + rdist2; % Space msd
         
                     % adding the term as part 2 need
                    % this saves all r^(t)-r0 at all times.
                    %part2rms(ip,jtM,1)=xy1;%-x0(ip); %HI
                    %part2rms(ip,jtM,2)=xy2;%-y0(ip); %HI
                    %part2rms(ip,jtM,3)=xy3;%-z0(ip); %HI
                    Wkinet(ip,jtM) = W;         % Evolution of kinetic energy
                   %is_in(ip,jtM) = 1;                                         %NS: the particle is inside
                    end
                end
               
            end
        end % EO monitoring
        
        %%
  %  end % EO in-the-box check
    end
   
    %% escape statistics
   
    if iBCs == 0 % Escape statistics; only for open BCs
        if lcollSwitch
            lFree = l * kFree;
            kstepColl(ip) = kstepColl(ip) + 1;
           
            N1 = randn;
            N2 = abs(randn) * sign(v);
            nuColl_v = nuColl .* vth.^3 ./ abs(v).^3;
            nuColl_arr(ip) = nuColl_arr(ip) + nuColl_v;
            rootCollFUN = @(tau) abs(StauFUN(tau, v, nuColl_v, N1, N2, kTm)) - lFree;
            tau1 = lFree / abs(v);
            tau2 = sqrt(lFree / abs(sqrt(kTm)*N1));
            if tau1 > tau2;
                tau12 = tau1;
                tau1 = tau2;
                tau2 = tau12;
            end
            iter1 = 0;
            iter2 = 0;
            root1 = rootCollFUN(tau1); % For fewer function calls in 2nd while
            while root1 > 0
                tau1 = 0.1 * tau1;
                iter1 = iter1 + 1;
                root1 = rootCollFUN(tau1);
            end
            root2 = rootCollFUN(tau2);
            while root2 < 0
                tau2 = 2 * tau2;
                iter2 = iter2 + 1;
                root2 = rootCollFUN(tau2);
            end
            tau = fzero(rootCollFUN,[tau1 tau2]);
            mu = exp(-nuColl_v * tau);
            vtau = v * mu + sqrt(kTm * (1-mu^2)) * N1;
           
            if lCstats
                iter1Arr = [iter1Arr; iter1];
                iter2Arr = [iter2Arr; iter2];
                tau1Arr = [tau1Arr; tau1];
                tau2Arr = [tau2Arr; tau2];
                root1Arr = [root1Arr; root1];
                root2Arr = [root2Arr; root2];
                tauArr = [tauArr; tau];
                rootsArr = [rootsArr; rootCollFUN(tau)];
                lFreeArr = [lFreeArr; lFree];
                vtauArr = [vtauArr; vtau];
                dvtauArr = [dvtauArr; vtau-v];
                nuCollArr = [nuCollArr; nuColl_v];
                muArr = [muArr; mu];
            end
            %v = abs(vtau);
            v = vtau;
            gamma = 1 ./ sqrt(1-(v./c).^2);
            if gamma > gammaLim
                W = (gamma-1) .* Wrest;
            else
                W = m*v.^2 / 2;
            end
           
            if gamma < 1    % Particle velocity is too low
                W = W_prev;
                nCignore(ip) = nCignore(ip) + 1;
                iCignore = [iCignore; ip];
                tCignore = [tCignore; t];
            end
            if abs(rootCollFUN(tau)) > 1.e-6
                rootInaccur = [rootInaccur; ip rootCollFUN(tau)];
            end
        else
            if t>= t_final
                             
                 t_esc(ip) = 0;
                 W_esc(ip) = 0;            % for the case where the particle is still inside the box
                
            else
                dt1  = tau_check(1,j-2);
            if x(1,j-2)<0
            
            dist(1,1)   = abs((dmin-x(1,j-2))/(sintheta(1,j-2)*cosphi(1,j-2)));
            speed(1,1)  = abs(x(1,j-1)-x(1,j-2))/dt1;
           else
             dist(1,1)  = abs((dmax-x(1,j-2))/(sintheta(1,j-2)*cosphi(1,j-2)));
             speed(1,1) = abs(x(1,j-1)-x(1,j-2))/dt1;
           end
         
         if y(1,j-2)<0
            dist(1,2)  = abs((dmin-y(1,j-2))/(sintheta(1,j-2)*sinphi(1,j-2)));
            speed(1,2) = abs(y(1,j-1)-y(1,j-2))/dt1;
         else
            dist(1,2)  = abs((dmax-y(1,j-2))/(sintheta(1,j-2)*sinphi(1,j-2)));
            speed(1,2) = abs(y(1,j-1)-y(1,j-2))/dt1;
         end
         
         if z(1,j-2)<0
            dist(1,3)  = abs((dmin_z-z(1,j-2))/(costheta(1,j-2)));
            speed(1,3) = abs(z(1,j-1)-z(1,j-2))/dt1;
         else
            dist(1,3)  = abs((dmax_z-z(1,j-2))/(costheta(1,j-2)));
            speed(1,3) = abs(z(1,j-1)-z(1,j-2))/dt1;
         end
       
         tau33=zeros(1,3);
         for k=1:3
             tau33(1,k)=abs(dist(1,k)/speed(1,k));  % Estimate the time between last
         end                                        % scattering event and boundary crossing
                                                    
         tau22=min(tau33);                          % Choose the minimum time and therefore the correct one
                                                       
           
             t_esc(ip) = t+tau22;                  % Add the estiomated time to the time up to that point
             W_esc(ip) = W;                        % Estimate escape energy
             %v_esc(ip) = v;
            end
        end % EO final collisions
      
    end
    if lscattSwitch
       if (kkicks(ip) ~= 0)
            mDelW_W(ip) = DelW_W ./ kkicks(ip); % < (\Delta W) / W > per particle
           % meanFP(ip) = l .* kFreePath(ip) ./ kkicks(ip);
        else
            mDelW_W(ip) = DelW_W;
           % meanFP(ip) = l .* kFreePath(ip);
        end
   end

    % Complete timeseries
    %if lTScomplete
        %W_Cell{ip} = W_Arr;
        dW_Cell{ip} = dW_Arr;
        %dW_Cell_overW(ip)= dW_Arr_overW;
       % t_Cell{ip} = t_Arr;
        %dtScat_Cell{ip} = dtScat_Arr;
        %kicks_Cell{ip} = kicks_Arr;
       % if imodel > 10, model_Cell{ip} = model_Arr; end
    %end
    
    %{
     if kkicks(ip)>=N  %re-run the code for this particle
        ip=ip-1;
        N_check=1;
    end
    %}
end % Particle loop
if lTime, elapsedTime = toc; end
%% ------------------------------------------------------------------------
%% ----- Saving data  % After this point it's just for saving the data
if lmonitSwitch
    % Space: mean displacements and diffusion coefficients
    %Drr_clas = 0.5*(rmsd1./is_in) ./tM;
   
   %save('r_coeffs006.mat',...
    %  'tM','ntM','tfin', 'kpM', 'rmsd', 'rmsd_', 'rmsd1', 'Drr_clas'...
     % );
   
    % Energies
    Wkinet_eV = Wkinet .* erg2eV;

  % save('W_evol006.mat',...
    %   'tM', 'Wkinet','erg2eV','Wkinet_eV'...
    % );
end

%% Escape data
if iBCs == 0 % Escape statistics; only for open BCs
    W_esc_eV = W_esc .* erg2eV;
  % save('esc006.mat',...
   %    't_esc', 'W_esc_eV', 'v_esc'...
   % );
end

if lscattSwitch
    %% Kicks etc.
    sumdWeV_p = sumdW_p .* erg2eV;
    
    W0_eV = W0 .* erg2eV;

clear cosphi sinphi sintheta costheta dr Wkinet W0 W_esc v0
clear cosphi costheta dr Wkinet W0 W_esc v0

disp('ok')

   % save('scattPathtperiod.mat',...
        % 'mDelW_W', 'kFreePath', 'meanFP',...
        % 'hdonkicks', 'forwkicks'; 'kkicks',...
         %'sumdWeV_p', 'countdW_p'...
        %);
    if lmonitSwitch
       % save('scattPathpaeriod.mat',...
             %'ktkicks',...
       % '-append');
    end
    %% Ignored scattering events
    %save('IgnoredScattperiod.mat',...
        % 'nSignore', 'iSignore', 'tSignore'...
       % );
end

%% Collision counters
if lcollSwitch
    nuCollMean = nuColl_arr ./ kstepColl;
    nuCollNum = sum(nuColl_arr) ./ sum(kstepColl);
    %save('collisionsperiod.mat',...
        % 'kstepColl', 'nuColl_arr', 'nuColl', 'nuCollNum', 'nuCollMean',...
        % 'nCignore', 'iCignore', 'tCignore', 'rootInaccur'...
       % ); % realy needed?
end
%% CS debug statistics
if lstats
    if imodel == 2
        %save('CSstattesetperiod.mat',...
           %  'nScat1', 'ngamma', 'maxgamma', 'gammaArr'...
           % );
        if lBdepend
          %  save('CSstatesettperiod.mat',...
             %    'BeffArr', 'EeffArr', 'LeffArr',...
           % '-append');
        else
            if lEplaw
              %  save('CSstatperiod.mat',...
                  %  'EeffArr',...
                    %'-append');
            end
            if lLeffrand
              %  save('CSstatperiod.mat',...
                %    'LeffArr',...
                 %   '-append');
            end
        end
        if idtCS == 1
          %  save('CSstatperiod.mat',...
         %        'ndtScat', 'maxdtScat', 'dtCSArr',...
         %   '-append');
        end
    end
end
%% Complete timeseries
if lTScomplete
   % save('WtTSeriesperiod.mat',...
         %'W_Cell', 'dW_Cell', 't_Cell', 'dtScat_Cell', 'kicks_Cell' );%', model_Cell' ...
        %);
end
%% Meta data
cwdFull = pwd;
[path2cwd, cwd, ~] = fileparts(cwdFull);
%save('metaperiod.mat',...
  %   'path2cwd', 'cwd'...
 %   );

 
if lTime
   % save('metaperiod.mat',...
       %   'elapsedTime'...
   % ,'-append');
end



