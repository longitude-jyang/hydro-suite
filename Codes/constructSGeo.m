function s_geo_v = constructSGeo(s_geo, s_mech,options)
% constructSGeo.m contruct the geometry as consist of strips
%     Input:
%      {s_geo} is a structure type array that contains all the required input 
%      for geometry buildings.
%         type   type of structure:  riser or floating wind turbine (FWT)
%         L1     length of structure above water level
%         L2     length of structure below water level
%         Li     length of strip (regular case)
%         D      diameters
%     Output:
%       {s_geo_v} that contains following fields/values, for all strips
%           [xi_v]   coordinate of strips
%           [Li_v]   length of strips
%           [Di_V]   outer diameter of strips
%           [Di_v_I] inner diameter of strips
%           [Mi_v]   structural mass of strips
%           [Mia_v]  added mass from fluid inertia of strips
%           (xc)     position of centre of gravity
%           (xcb)    position of centre of buoyancy
%           (xm)     position of mooring line spring (if any)

% 07/07/2020 @ Franklin Court, Cambridge  [J Yang] 
% 09/04/2021 @ Franklin Court,, Cambridge  [J Yang]  --> code checked
% 18/01/2022 @ Franklin Court,, Cambridge  [J Yang]  --> add output for
% Bi_v (for Ca=1, Mia_v = Bi_v, but in random version, it is different)

    isModelTest = options.isModelTest;

    L1      = s_geo.L1;

    rho     = s_mech.rho(1);
    rho_f   = s_mech.rho(2);
    rho_oil = s_mech.rho(3);

    Ca      = s_mech.Morrison(1);

    Nsec    = s_geo.Nsec;  % number of sections
    N       = s_geo.N; % total number of strips

    % initialise the vectors 
    
    % vector for mass of all strips
    Mi_v    = zeros(N,1); 
    % vector for added mass (the ones above fluid surface will be zeros
    Mia_v   = zeros(N,1); 
    % vector for buyancy
    Bi_v    = zeros(N,1);
    % vector for 2nd moment of area
    Ii_v    = zeros(N,1); 
    
    % vector for coordinate
    xi_v    = zeros(N,1);
    % vector for lengths
    Li_v    = zeros(N,1);
    
    % vector for diameter
    Di_v    = zeros(N,1);
    Di_v_I  = zeros(N,1);
    Di_v_out = zeros(N,1);

    % for each section 
    for ss = 1 : Nsec

        Li = s_geo.Li(ss);                          
        Ns = s_geo.Ns(ss);

        Ls = s_geo.L(ss);
        
        rho=rho*s_geo.rhoC(ss);

        
        D    = s_geo.D (: , ss);
        Din  = D(1);  % inner diameter
        Dout = D(2);  % outer diameter
        De   = D(3);  % hydrodynamic equivalent diameter
        

        if ss == 1
            xi_v_s = linspace(0+Li/2,0+Ls-Li/2,Ns).';
        else
            xi_v_s = linspace(0+Li/2,Ls-Li/2,Ns).'...
                     + sum(s_geo.L(1:ss-1));
        end

        % for each strip 
        for ii=s_geo.Nsv{ss}

            % diameter of ith strip
            Di_v(ii)     = De; 
            Di_v_out(ii) = Dout; 
            Di_v_I(ii)   = Din; 

            % volume of ith strip
            Vout         = pi*(Dout/2)^2*Li; 
            Vin          = pi*(Din/2)^2*Li;
            Ve           = pi*(De/2)^2*Li;

            %  mass of ith strip (structure + oil)
            Mi           = rho*(Vout-Vin)+rho_oil*Vin ; 

            % 2nd moment of area, for circular cross-section pi*radius^4/4
            Ii           = pi*((Dout/2)^4-(Din/2)^4)/4; 

            
            % buoyancy of ith trip ( equal to the mass of displaced fluid)
            Bi           = Ve*rho_f; 
            % added mass from fluid, Bi is buoyancy of ith trip
            Mia          = Ca*Bi; 

            % put them into vector form
            Mi_v(ii)     = Mi;
            Mia_v(ii)    = Mia;
            Bi_v(ii)     = Bi;

            Ii_v(ii)     = Ii;
            Li_v(ii)     = Li;

        end

        xi_v(s_geo.Nsv{ss}') = xi_v_s;
    end

    % water level is zero for x coordinate
    xi_v            = xi_v-L1;
    
    % no fluid above fluid surface
    Mia_v(xi_v<0)   = 0;   
    Bi_v(xi_v<0)    = 0;  

    if isModelTest == 1 
        mb   = s_geo.mb;
        L_b  = s_geo.L_b;
        d    = s_geo.d;
        L_S  = s_geo.L(2);
        
        [~, mb_index] = min(abs(xi_v - ((d - L_S) - L_b)));
        
        Mi_v (mb_index) = Mi_v(mb_index) + mb;
        
    end
    
    % get centre of gravity , centre of buoyancy
    xc              = sum(Mi_v.*xi_v)/sum(Mi_v); 
    xcb             = sum(Mia_v.*xi_v)/sum(Mia_v); 

    % coor for spring (dummy spring?) 
    xm              = xi_v(end-3);          
    
    
    %-------------------------------------------------------------------------
    % put everything into structure array for output
    %-------------------------------------------------------------------------

    s_geo_v.N       = N;
    s_geo_v.x       = xi_v;
    s_geo_v.L       = Li_v;
    s_geo_v.D       = [Di_v_I Di_v_out Di_v];
    s_geo_v.m       = Mi_v;
    s_geo_v.ma      = Mia_v;
    s_geo_v.Bi      = Bi_v;
    s_geo_v.xc      = xc;
    s_geo_v.xcb     = xcb;
    s_geo_v.xm      = xm;

    s_geo_v.I       = Ii_v;

end
