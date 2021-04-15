function options = initialiseOpts (ModPar)
% initialiseOpts.m initialise the option parameters

% 07/07/2020 @ Franklin Court, Cambridge  [J Yang] 
% 09/04/2021 @ Franklin Court,, Cambridge  [J Yang]  --> code checked

    options.wavetype    = ModPar.wavetype;  % wave type: harmonic or random
    options.wavecode    = ModPar.wavecode;  % type of random waves
    options.Struct_type = ModPar.Struct_type;  % FWT or Riser
    options.bodyflex    = ModPar.bodyflex;  % body flexibility, rigid for FWT, flexible for riser
    options.TF          = ModPar.calcType;  %  time or frequency analysis
    options.Ndof        = ModPar.Ndof;      % number of shape functions
    options.Nvb         = ModPar.Nvb;       % number of modes to be kept
    options.Nfft        = ModPar.Nfft;      % N for fft
    options.Nwin        = ModPar.Nwin;      % number of windows for random TD
    options.wLim        = ModPar.wLim;      % Jonswap related parameter
    options.typeIntTLP  = ModPar.typeIntTLP;
    options.g           = ModPar.g; 
    
    options.isModelTest = ModPar.isModelTest;  % is this model test (wave tank)
    options.om_range    = ModPar.om_range;     % frequency range for harmonic analysis
    options.tmax        = ModPar.tmax;         % tmax for time analysis
    options.dtScale      = ModPar.dtScale;     %  scale resolution for time data

