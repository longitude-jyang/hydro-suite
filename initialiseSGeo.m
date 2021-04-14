function s_geo = initialiseSGeo (ModPar)
% initialiseSGeo.m initialise the geometric parameters

% 07/07/2020 @ Franklin Court, Cambridge  [J Yang] 
% 09/04/2021 @ Franklin Court,, Cambridge  [J Yang]  --> code checked

    % L1: length above water; L2: length below water; L: total length
    s_geo.L1   = ModPar.L1;     
    s_geo.L2   = ModPar.L2;
    
    % if model testing
    s_geo.d    = ModPar.d;   % d: wave tank depth
    s_geo.L_b  = ModPar.L_b; % ballast position
    s_geo.mb   = ModPar.mb;  % ballast weight
    
    % number of sections, each section could have different geometries 
    s_geo.Nsec = ModPar.Nsec; % number of sections
    
    % strip parameters, should be vectors with size corresponding to the 
    % number of sections. Li: length of each strip;  Ns: number of strips 
    % of each section; N: total number for all sections
    s_geo.L    = ModPar.L;
    s_geo.Li   = ModPar.Li;
    s_geo.Ns   = s_geo.L./s_geo.Li; 
    s_geo.N    = sum(s_geo.Ns);
    
    % put index of strips into a vector, corresponding to each section 
    s_geo.Nsv  = cell (s_geo.Nsec , 1);
    s_geo.Nsv{1} = 1 : s_geo.Ns(1);
    
    if  s_geo.Nsec >= 2  
        for iSec = 1 : s_geo.Nsec - 1  
            s_geo.Nsv { iSec + 1 } = sum (s_geo.Ns(1 : iSec)) + 1 ...
                                 : sum (s_geo.Ns(1 : iSec + 1));
        end
    end
  
    % diameter, vector with size of Nsec 
    Din  = ModPar.Din;
    Dout = ModPar.Dout;
    De   = ModPar.De;
    
    % density coefficient for each section, size of Nsec
    s_geo.rhoC = ModPar.rhoC;
    
    % each column correspond diameters (inner, outer, equivalent) of each
    % section
    s_geo.D = [Din ; Dout ; De];     