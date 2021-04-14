function s_mech = initialiseSMech ( ModPar )
% initialiseSMech.m  initialise the mechanical parameters 

% 07/07/2020 @ Franklin Court, Cambridge  [J Yang] 
% 09/04/2021 @ Franklin Court,, Cambridge  [J Yang]  --> code checked
    
    % Morison parameters, Ca/Cm:added mass coef; Cd:drag coef;
    Ca = ModPar.Ca; 
    Cd = ModPar.Cd;
    Cm = Ca+1; 
    
    % Material, E: Young's modulus; rho: material density
    rho      = ModPar.rho;
    rho_f    = ModPar.rho_f;
    rho_oil  = ModPar.rho_oil;
    E        = ModPar.E;

    % Top platform (if any) parameters, T0: top tension; TLP: platform
    % frequency dependent motion parameters
    T0        = ModPar.T0;
    T_TLP     = ModPar.T_TLP;
    H_TLP     = ModPar.H_TLP;
    gamma_TLP = ModPar.gamma_TLP;
    
    s_alpha   = ModPar.s_alpha;
    
    % for harmonic wave
    aw        = ModPar.aw ; 
    
    % parameters for random wave
    % significant wave height [meter] % wave periods [seconds]
    JONSWAP_w.Hs = ModPar.jonswapHs; 
    if isfield(ModPar, 'jonswapTz')
        JONSWAP_w.Tz = ModPar.jonswapTz; 
    end
    
    
    %-------------------------------------------------------------------------
    % put mechanical properties into structure array for output
    %-------------------------------------------------------------------------
    s_mech.E        = E;
    s_mech.rho      = [rho rho_f rho_oil];
    s_mech.Morrison = [Ca Cm Cd];
    s_mech.T0       = T0;
    s_mech.waveR    = JONSWAP_w;                % for random wave
    s_mech.waveH    = aw;                       % for harmonic wave
    s_mech.spring   = s_alpha;

    riser.H_TLP     = H_TLP;
    riser.T_TLP     = T_TLP;
    riser.gamma_TLP = gamma_TLP;

    s_mech.riser    = riser;