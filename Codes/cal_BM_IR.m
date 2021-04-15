function BM_IR = cal_BM_IR(s_geo_v, s_mech, wave, theta, theta2, response, options)
% cal_BM_IR calculates the bending moment with inertia relief (IR)
% this is mainly used for rigid bodys, where the 2nd derivatives are not
% available 

% 15/03/2019 @ JDG-01, Cambridge  [J Yang] 
% 09/04/2021 @ Franklin Court,, Cambridge  [J Yang]  --> code checked

    g = options.g;
    wavetype = options.wavetype;
    TF = options.TF;

    % geometry along the structure 
    N    = s_geo_v.N;
    xi_v = s_geo_v.x;
    Li_v = s_geo_v.L;
    Di_v = s_geo_v.D(:,3);
    
    % distribution of mass along the structure 
    Mi_v  = s_geo_v.m;
    Mia_v = s_geo_v.ma;

    % Morrison's parameters 
    Cm  = s_mech.Morrison(2);
    Cd  = s_mech.Morrison(3);
    
    % mooring spring 
    xm      = s_geo_v.xm;
    s_alpha = s_mech.spring;

    % top tension 
    T0 = s_mech.T0;

    % water density 
    rho_f = s_mech.rho(2);

    % wave amplitdue 
    aw   =  s_mech.waveH;
    
    
% -------------------------------------------------------------------
% first extract wave and structure responses 

    if TF == 1
        
        y  = response.frf;
        om_range  = wave.F_f;

        Nom = length(om_range);
        Mb  = zeros(N,Nom); % predefine empty matrix here

        % wave velocity        
        k_w = om_range.^2/g;  % wavenumber (based on linear wave theory)

        u_v     = 1i*om_range.*aw.*exp(-xi_v*k_w);  % fluid velocity, matrix format, N x Nom
        u_v_dot = -om_range.^2.*aw.*exp(-xi_v*k_w); % fluid acceleration,  N x Nom

        u_v(xi_v<0,:) = 0;
        u_v_dot(xi_v<0,:) = 0;

        % structural velocity
        u_s      = y.*(1i*om_range); % structure velocity
        u_s_dot  = -y.*om_range.^2; % structure acceleration

    elseif TF == 2

        y = response.yf;

        tf  = response.tf;
        Ntf = length(tf);
        Mb  = zeros(N,Ntf);% predefine empty matrix here

        % wave velocity
        u_v      = response.u;
        u_v_dot  = response.u_d;
        
        % structural velocity
        u_s     = response.yf_1d;
        u_s_dot = response.yf_2d;

    end
    
% -------------------------------------------------------------------
% set up relative velocity and linearise if TF == 1 (frequency domain)

    if TF==1
       
        % linearization coefficient
        if wavetype == 1
            
            gamma = 8/(3*pi)*abs(u_v-u_s); % gamma is now in matrix form
            
        elseif wavetype == 2
            
            dom = om_range(2) - om_range(1); % delta omega

            Sxx = wave.Sxx.';
            Sur = abs(u_v-u_s).^2.*Sxx;
            gamma = sqrt(8/pi)*sqrt(sum(Sur*dom,2)); % linearised coefficient, initialise (sqrt(8/pi)*std(ur))

        end

    elseif TF==2
        
        gamma = abs(u_v-u_s);  % call it gamma for convenience, but no linearization for time domain analysis
    
    end

% -------------------------------------------------------------------
% forces 

    % wave force 
    Fwave = (pi/4*(Di_v).^2*rho_f*Cm.*Li_v).*u_v_dot +...
        1/2*rho_f*Cd.*Di_v.*Li_v.*gamma.*((u_v-u_s)); %  N x Nom
    Fwave(xi_v<0,:) = 0; % only effective below fluid surface

    % inertial force
    Finer = (Mi_v + Mia_v).*u_s_dot; % F=ma

    % forces due to effective tension
    dTe = (Mi_v - Mia_v)*g;
    F_Te = (T0 - sum(triu(repmat(dTe,1,N))).').*theta2.*Li_v;

    % effective weight
    F_we = (dTe./Li_v).*theta.*Li_v;

    % spring force
    Kspring = zeros(N,1);
    Kspring(xi_v==xm) = s_alpha;
    Fspring = Kspring.*y;

    % overall force at each strip
    Fi_v = Fwave - Finer - Fspring + F_Te - F_we;

    % -------------------------------------------------------------------
        % boundary forces, as external force, needs to be considered for fixed/pinned ends 
        % with free boundary, the boundary forces should be small (in terms of overall contribution to moment)
        
        Fb1 = sum(Fi_v.*(xi_v(end) - xi_v),1)/(xi_v(end)-xi_v(1)); % force at the boundary (top)
        Fb2 = sum(Fi_v.*(xi_v(end) - flipud(xi_v)),1)/(xi_v(end)-xi_v(1)); % force at the boundary (bottom)
       
        Fi_v(1,:)   = Fi_v(1,:)   - Fb1;
        Fi_v(end,:) = Fi_v(end,:) - Fb2;
        
    % -------------------------------------------------------------------
    % Bending Moment

    for ii=1:N
        Mb(ii,:) = sum(Fi_v(1:ii,:).*(xi_v(1:ii) - xi_v(ii)),1);
    end

    % -------------------------------------------------------------------
    % Output results

    BM_IR.Smm = Mb;

end

    
