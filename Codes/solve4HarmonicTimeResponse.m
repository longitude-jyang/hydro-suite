function response = solve4HarmonicTimeResponse( s_geo_v, s_mech, wave, mk, mode, options)
% solve4HarmonicTimeResponse.m solve the structural response due to
% harmonic wave exctations, in time domain 
% note that the variable at Line78: number of cycles to be used 

% 23/07/2020 @ Franklin Court, Cambridge  [J Yang] 
% 09/04/2021 @ Franklin Court,, Cambridge  [J Yang]  --> code checked

    % initialise inputs 
    Struct_type = options.Struct_type;
    g           = options.g;
    
    N           = s_geo_v.N;
    xi_v        = s_geo_v.x;
    aw          = s_mech.waveH;
    
    Phi_x       = mk.Phi_x;
    Phi_x_0d    = Phi_x.d0;
    
    Nvb         = mode.Nvb;
    Vb          = mode.Vb;
    
    Nt          = wave.Nt;
    tspan       = wave.t;
    tmax        = wave.tmax;

    om_range_t  = wave.T_f;
    N_om_t      = length(om_range_t);
                
    TLP_H       = wave.TLP;
    
    % initialise empty matrices 
    y_peak      = zeros(N,N_om_t);           % peak displacement 
    
    SmmEI       = zeros(N,N_om_t);           % bending moment (EI method)
    SigmaEI     = zeros(N,N_om_t);           % bending stress (EI method)

    SmmIR       = zeros(N,N_om_t);           % (IR method)
    SigmaIR     = zeros(N,N_om_t);
                
    %  loop over frequency range 
    jj  =  0;
    for om = om_range_t
        jj = jj+1;
        
        % prescribed top motion, its 1st and 2nd derivatives
        q_1     = TLP_H.q_1(jj,:);
        q_1_1d  = TLP_H.q_1_1d(jj,:);
        q_1_2d  = TLP_H.q_1_2d(jj,:);
                    
        wave.q_1_om = [q_1;q_1_1d;q_1_2d];
        
        %---------------------------------------
        % ODE stepping for time domain solutions
        phi   = 0; % initial phase
        a0    = zeros(Nvb*2,1); % initial value for a 
        [t,a] = ode45(@(t,a) odeFloatCylinder_beam_v2(t,a,s_geo_v,s_mech,mk,mode,wave,options,phi,om,aw,Phi_x_0d), tspan, a0);
                    
        %---------------------------------------            
        % mode acceleration
        [q,a_2d] = mode_acc(s_geo_v,s_mech,mk,mode,wave,options,a,t,om,aw,phi,Phi_x_0d);      
        
        %---------------------------------------
        % transform response from generalized coordinates to physical coordinates 
        if Struct_type == 1
            y    = Phi_x_0d*q; 
            y_1d = Phi_x_0d*Vb*a(:,Nvb+1:end).';
            y_2d = Phi_x_0d*Vb*a_2d; % acceleration

        elseif Struct_type == 2
            y    = Phi_x_0d(:,2:end)*q+Phi_x_0d(:,1)*q_1;
            y_1d = Phi_x_0d(:,2:end)*Vb*a(:,Nvb+1:end).'+Phi_x_0d(:,1)*q_1_1d;
            y_2d = Phi_x_0d(:,2:end)*Vb*a_2d+Phi_x_0d(:,1)*q_1_2d; 
        end
        
        %-----------------------------------------------------------
        % for harmonic excitation, instead of fft, it is much easier to cut the response , take a few cycles and get peak2peak value
        [~,itm]      = min(abs(tspan-(tmax-2*pi/(om)*5)));  % take 5 cycles <---------------
        tcut         = itm : Nt;
        tf           = tspan(tcut);
        yf           = y(:,tcut);
        y_peak(:,jj) = max(yf,[],2); % get peak value for all columns, representing frequency amplitude

        if Struct_type == 1
            qf = q(:,tcut);
        elseif Struct_type == 2
            qf = [q_1(tcut);q(:,tcut)];
        end
        
        yf_1d = y_1d(:,tcut);
        yf_2d = y_2d(:,tcut);

        % -------------------------------------------------------------- 
        % before saving, also compute particle velocity and acceleration 
        [u_v,u_v_dot] = genTimeHarmonicWave (aw,om,tspan,phi,xi_v,g); 

        
        %-------------------------------------------------
        % put the responses in structure array 
        response.tf     = tf;              % time vector
        response.yf     = yf;              % displacement vector
        response.yf_1d  = yf_1d;              % 1st derivative 
        response.yf_2d  = yf_2d;              % 2nd derivative
        response.qf     = qf;              % generalised response 

        response.u      = u_v(:,tcut);
        response.u_d    = u_v_dot(:,tcut);

        %-------------------------------------------------
        % estimate bending moment and stress 
        BM  = cal_BM_AC(s_geo_v, s_mech, mk, wave, response, options);

        SmmEI(:,jj)   = BM.SmmEI;
        SigmaEI(:,jj) = BM.SigmaEI;

        SmmIR(:,jj)   = BM.SmmIR;
        SigmaIR(:,jj) = BM.SigmaIR;
    end

    response.yt     =  y;
    response.y_peak =  y_peak;

    BM.SmmEI   = SmmEI;
    BM.SigmaEI = SigmaEI;

    BM.SmmIR   = SmmIR;
    BM.SigmaIR = SigmaIR;

    response.BM = BM;