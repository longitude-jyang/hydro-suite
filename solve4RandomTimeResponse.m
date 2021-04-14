function response = solve4RandomTimeResponse( s_geo_v, s_mech, wave, mk, mode, options)
% solve4RandomTimeResponse.m solve the structural response due to
% random wave exctations, in time domain 

% 23/07/2020 @ Franklin Court, Cambridge  [J Yang] 

    Struct_type = options.Struct_type;
    g           = options.g;
    Nwin        = options.Nwin;

    Nt       = wave.Nt;
    tspan    = wave.t;
    om_range = wave.F_f;
    
    Nvb      = mode.Nvb;
    Vb       = mode.Vb;
    
    Phi_x    = mk.Phi_x;
    Phi_x_0d = Phi_x.d0;
                
    TLP_R    = wave.TLP;
    q_1      = TLP_R.q_1.';
    q_1_1d   = TLP_R.q_1_1d.';
    q_1_2d   = TLP_R.q_1_2d.';
    
    aw       = s_mech.waveH;
    

    
    %---------------------------------------
    % ODE stepping for time domain solutions
    om  = []; % not used
    phi = [];
    a0       = zeros(Nvb*2,1); % initial value
    [t,a]    = ode45(@(t,a) odeFloatCylinder_beam_v2(t,a,s_geo_v,s_mech,mk,mode,wave,options,phi,om,aw,Phi_x_0d), tspan, a0);
 
    %---------------------------------------            
    % mode acceleration
    [q,a_2d] = mode_acc(s_geo_v,s_mech,mk,mode,wave,options,a,t,om,aw,phi,Phi_x_0d);
    
    %---------------------------------------
    % transform response from generalized coordinates to physical
    % coordinates 
    if Struct_type == 1
        y    = Phi_x_0d*q; 
        y_1d = Phi_x_0d*Vb*a(:,Nvb+1:end).';
        y_2d = Phi_x_0d*Vb*a_2d; % acceleration

    elseif Struct_type == 2
        y    = Phi_x_0d(:,2:end)*q+Phi_x_0d(:,1)*q_1;
        y_1d = Phi_x_0d(:,2:end)*Vb*a(:,Nvb+1:end).'+Phi_x_0d(:,1)*q_1_1d;
        y_2d = Phi_x_0d(:,2:end)*Vb*a_2d+Phi_x_0d(:,1)*q_1_2d; 
    end

    %-----------------------------------
    % cut response to remove transients
    tcut = Nt-round(Nt/1.2):Nt;

    yf    = y(:,tcut).';
    yf_1d = y_1d(:,tcut);
    yf_2d = y_2d(:,tcut);

    if Struct_type == 1
        Hq = q;
    elseif Struct_type == 2
        Hq = [q_1; q];
    end
                
    qf  = Hq(:,tcut);
    tf  = t(tcut);
    NtN = length(tf);
         
    %----------------------------------
    % convert to PSD
    fs       = om_range(end)/(2*pi)*2;
    [Ayy,f]  = pwelch(yf,round(NtN/Nwin),[],[],fs) ;
    Syy      = Ayy.'/(2*pi);

    rms_Syy  = sqrt(sum(Syy*2*pi*(f(2)-f(1)),2)); % std value for all strips
                
    % save for  output
    response.fSyy    = f;
    response.Syy     = Syy;
    response.rms_Syy = rms_Syy;
    response.rms_yf  = std(yf);

    response.tf      = tf;
    response.yf      = yf.';
    response.yf_1d   = yf_1d;
    response.yf_2d   = yf_2d;

    response.qf      = qf;

    % also save cutted wave velocity for BM recover
    u_v              = wave.u;
    u_v_dot          = wave.u_d;

    response.u       = u_v(tcut,:).';
    response.u_d     = u_v_dot(tcut,:).'; 