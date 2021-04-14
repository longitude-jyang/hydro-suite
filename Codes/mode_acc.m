function [q_v,a_2d] = mode_acc(s_geo_v,s_mech,mk,mode,wave,options,a,t,om,aw,phi,Phi_x)
% mode_acc recovers the response with mode acceleration approach 
% this code is quite similar to odeFloatCylinder_beam_v2.m 
% main differences are:
%  -  no time interpolation here because time vectors are as defined
%  -  has to loop through time vector 

% 15/03/2019 @ JDG-01, Cambridge  [J Yang] 
% 09/04/2021 @ Franklin Court,, Cambridge  [J Yang]  --> code checked

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
    M_1   = mk.constrainM;
    K_1   = mk.constrainK;
    Ndof  = mk.Ndof;
    
    % modes
    Vb     = mode.Vb;   
    Nvb    = mode.Nvb;
    D_om_b = mode.omb;
    

    % Morrison's coefficients 
    Cm  = s_mech.Morrison(2);
    Cd  = s_mech.Morrison(3);
        
    t = t'; % put t as row vector   
    if wavetype == 1  % harmonic design wave 
        
            [u_v,u_v_dot] = genTimeHarmonicWave (aw,om,t,phi,xi_v,g); 
            
            if Struct_type == 2 % for riser, there is presecried top motion

                q_1    = wave.q_1_om(1,:);
                q_1_1d = wave.q_1_om(2,:);
                q_1_2d = wave.q_1_om(3,:);
            end
          
    elseif wavetype == 2 % random wave

            u_v      = wave.u.';% convert to column vector, each colum represent the structure at one of t , N x Nt
            u_v_dot  = wave.u_d.';
            
            if Struct_type == 2 %
                TLP_R = wave.TLP;
                
                q_1    = TLP_R.q_1.';
                q_1_1d = TLP_R.q_1_1d.';
                q_1_2d = TLP_R.q_1_2d.';
            end
    end 
    
    % ------------------------------------
    % response calculated from ode45  
    a_0d = a(:,1:Nvb).';
    a_1d = a(:,Nvb+1:end).';
     
    if  Struct_type == 1    
        
        ur_v = u_v - Phi_x*Vb*a_1d; % relative velocity
                
    elseif Struct_type == 2
        
        Phi_x_1 = Phi_x(:,1);    
        ur_v    = u_v - Phi_x(:,2:end)*Vb*a_1d - Phi_x_1*q_1_1d; % relative velocity
    end
    
    
    % run loop through the time vector from here (can only do matrix calculation at each t)
    Nt   = length(t);
    q_v  = zeros(Ndof,Nt);
    a_2d = zeros(Nvb,Nt);
    
    for ii = 1 : Nt
        
        % damping matrix  
        BDi_v = 1/2*rho_f*Di_v.*Li_v.*Cd.*abs(ur_v(:,ii)); 
        BDi_v(xi_v<0) = 0; % only effective below fluid surface 

        if Cd == 0  % if no drag force, add artificial damping 
             C = 0.01*M;
             disp('If Cd = 0, artificial damping added')
        else
             % Mathematically we form C as : C=(Phi_x.')*diag(BDi_v)*Phi_x, but
             % since the middle one is diagonal, we can use the below vector
             % form to same some space ( in case the diagonal matrix is big)
             C = (Phi_x.').*BDi_v.'*Phi_x; 
        end


        % force vector 
        FMi_v = pi/4*(Di_v).^2*rho_f*Cm.*Li_v.*u_v_dot(:,ii) +...
             1/2*rho_f*Cd.*Di_v.*Li_v.*u_v(:,ii).*abs(ur_v(:,ii)); 
        FMi_v(xi_v<0) = 0; % only effective below fluid surface 

        F = (Phi_x.')*FMi_v;


        if Struct_type == 2 % for riser, there is presecried top motion 

            C_1 = C(2:end,1);
            C   = C(2:end,2:end);
            F   = F(2:end);

            F   = F - M_1*q_1_2d(:,ii) - K_1*q_1(:,ii) - C_1*q_1_1d(:,ii);
        end

        %----------------------------------------
        % calculate response at in modal coodinates 
        Fb = Vb.'*F;  % force in modal coordinates 
        Cb = Vb.'*C*Vb; % damping in modal coordinates 

        a_2d(:,ii) = Fb - Cb*a_1d(:,ii) - diag(D_om_b)*a_0d(:,ii);

        % mode acceleration 
        q_v(:,ii) = K\(F - M*Vb*a_2d(:,ii) - C*Vb*a_1d(:,ii)); % response in generalised coodindates 
    end