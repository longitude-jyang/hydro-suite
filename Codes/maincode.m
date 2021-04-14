function Out = maincode (ModPar)
% maincode.m is the main program for the suite of code to calculate the
% hydrodynamic response of a simple floating structure 
% refer to code manual for the details of input and output 

% 07/07/2020 @ Franklin Court, Cambridge  [J Yang] --> sort out previous
% versions into this format

% 07/04/2021 @ Franklin Court, Cambridge  [J Yang] --> run checks with example 
% caluclation for both FWT and riser 


%-----------------------------------------
% Section A: initialise 
%-----------------------------------------

    % (A-1)  initialise mechanical parameters  
           s_mech = initialiseSMech (ModPar) ;

    % (A-2)  initialise geometrical parameters
           s_geo  = initialiseSGeo (ModPar) ; 

    % (A-3)  initialise options 
           options = initialiseOpts (ModPar) ;
           
%-----------------------------------------
% Section B: construct       
%-----------------------------------------

    % (B-1)  construct geometries of the strutures 
           s_geo_v    = constructSGeo(s_geo, s_mech,options) ;     
           
    % (B-2)  construct M-K matrix
           [mk, mode] = constructMK(s_geo, s_geo_v, s_mech, options);
            
    % (B-3)  construct ocean waves 
            wave      = constructWave (s_geo_v, s_mech, options);

%-----------------------------------------
% Section C: solve       
%----------------------------------------- 

    % (C-1) Solve in time or frequency domain  
    response = solve4Response(s_geo, s_geo_v, s_mech, wave, mk, mode, options);

    % For harmonic wave in time domain, the calculation needs to be done at
    % each frequency and therefore it was done in the previous response computation
    % at the same time; for other types, here we extract the bending
    % moment/stress 
    if options.TF == 2 && options.wavetype == 1  
        BM = response.BM;
    else
        BM = cal_BM_AC(s_geo_v, s_mech, mk, wave, response, options);
    end

    
%-----------------------------------------
% Section D: assemble       
%-----------------------------------------    
    Out.BM = BM;
    Out.response = response; 

    switch ModPar.wavetype
        case 2
            Out.rms_Syy = response.rms_Syy;
    end

    % output frequency vector 
    switch ModPar.calcType
        case 1
            Out.freq = wave.F_f;
        case 2
            switch ModPar.wavetype
                case 1
                    Out.freq = wave.T_f;
                case 2
                    % Why 2*pi?
                    Out.freq = response.fSyy.'*2*pi;
            end
    end

    % output coordinate 
    Out.xCoord = s_geo_v.x;
    Out.wave = wave;
    Out.mode = mode;

end

