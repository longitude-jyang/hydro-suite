% odeFloatCylinder_beam solves floatcylinder_beam to a harmonic wave using ode
% function in the time domian; 
% the equation of motion if coverted to 1st order
% the drag force from Morrison's equation is kept as nonlinear 
function a_dot=odeFloatCylinder_beam_v2(t,a,s_geo_v,s_mech,mk,mode,wave,options,phi,om,aw,Phi_x)

% 15/03/2019 @ JDG-01, Cambridge  [J Yang] 
% 01/07/2019 @ JDG-01, Cambridge  [J Yang] --> update to compute the
% response with modal truncation 
% 07/04/2021 @ Franklin Court, Cambridge  [J Yang] --> code checked 

    g = options.g;
    wavetype = options.wavetype;
    Struct_type = options.Struct_type;

    % geometry 
    xi_v = s_geo_v.x;
    Li_v = s_geo_v.L;
    Di_v = s_geo_v.D(:,3);
    rho_f = s_mech.rho(2);
    
    % mass and stiffness matrices
    M = mk.M;
    K = mk.K;
    M_1 = mk.constrainM;
    K_1 = mk.constrainK;
    
    % modes
    Vb     = mode.Vb;   
    Nvb    = mode.Nvb;
    D_om_b = mode.omb;

    % Morrison's coefficients 
    Cm  = s_mech.Morrison(2);
    Cd  = s_mech.Morrison(3);
    

    if wavetype == 1  % harmonic design wave 
        
            [u_v,u_v_dot] = genTimeHarmonicWave (aw,om,t,phi,xi_v,g); 
            
            % for riser, there is presecried top motion   
            if Struct_type == 2  
                
                tspan = wave.t;
                q_1   = interp1(tspan,wave.q_1_om(1,:),t);
                q_1_1d = interp1(tspan,wave.q_1_om(2,:),t);
                q_1_2d = interp1(tspan,wave.q_1_om(3,:),t);
                
            end
          
    elseif wavetype == 2 % random wave

            % pick up pre-generated random wave time history 
            u_v_1 = wave.u; 
            u_v_dot_1 = wave.u_d;
            tspan = wave.t;
                   
            % interp to get the correponding excitation at time t, u_v is now a row vector representing response of all xi @ t
            u_v     = interp1(tspan,u_v_1,t);  
            u_v_dot = interp1(tspan,u_v_dot_1,t);
            u_v     = u_v.'; % convert to column vector 
            u_v_dot = u_v_dot.';
            
            % above water, switch off fluid data
            u_v(xi_v<0)     = 0; 
            u_v_dot(xi_v<0) = 0; 
            
            % for riser, there is presecried top motion   
            if Struct_type == 2 
                TLP_R = wave.TLP;
                
                q_1    = interp1(tspan,TLP_R.q_1,t);
                q_1_1d = interp1(tspan,TLP_R.q_1_1d,t);
                q_1_2d = interp1(tspan,TLP_R.q_1_2d,t);
            end
    end 
    
    % define relative velocity between the structure and the fluid 
    if  Struct_type == 1    
        
        ur_v = u_v - Phi_x*Vb*a(Nvb+1:end);   
        
    elseif Struct_type == 2
        
        Phi_x_1 = Phi_x(:,1); % for riser, there is presecried top motion             
        ur_v    = u_v - Phi_x(:,2:end)*Vb*a(Nvb+1:end) - Phi_x_1*q_1_1d;
    end
    

    % form damping matrix  (it is dependent on relative velocity as from Morrison'equation)
    BDi_v = 1/2*rho_f*Di_v.*Li_v.*Cd.*abs(ur_v); % damping vector 
    BDi_v(xi_v<0) = 0; % only effective below fluid surface 

    if Cd == 0  % if no drag force, add artificial damping 
         C = 0.01*M;
    else         
         % Mathematically we form C as : C=(Phi_x.')*diag(BDi_v)*Phi_x, but
         % since the middle one is diagonal, we can use the below vector
         % form to same some space ( in case the diagonal matrix is big)
         C = (Phi_x.').*BDi_v.'*Phi_x; 
    end
    

    % form force vector   
    FMi_v = pi/4*(Di_v).^2*rho_f*Cm.*Li_v.*u_v_dot +...  % inertia part
         1/2*rho_f*Cd.*Di_v.*Li_v.*u_v.*abs(ur_v);       % drag part, dependent on relative velocity as from Morrison'equation
    FMi_v(xi_v<0) = 0; % only effective below fluid surface 

    F = (Phi_x.')*FMi_v; % transform into shape function coordinate 

    
    if Struct_type == 2 % for riser, there is presecried top motion 

        C_1 = C(2:end,1);
        C   = C(2:end,2:end);
        F   = F(2:end);

        F   = F - M_1*q_1_2d - K_1*q_1 - C_1*q_1_1d;
    end
    
    
    %----------------------------------------
    % in modal coodinates (mass matrix is identity matrix, stiffness matrix is diagonal with eigenvalues on the diagonal)
    Fb = Vb.'*F;    % force in modal coordinates 
    Cb = Vb.'*C*Vb; % damping in modal coordinates 
 
    % populate A and B matrix for 1st order equation of motion and prepare for ode time integrations 
    Bp = [zeros(Nvb,Nvb); eye(Nvb,Nvb)]*Fb;

    A  = [zeros(Nvb,Nvb) -eye(Nvb,Nvb);...
          diag(D_om_b) Cb];

    a_dot = Bp - A*a; % 1st derivative equation to solve 
