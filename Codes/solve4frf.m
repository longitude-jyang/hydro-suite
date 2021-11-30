function [H, Hq] = solve4frf(s_geo_v, s_mech, mk, mode, wave, options)
% solve4frf.m solve frequency response fruntion (frf) 

% 23/07/2020 @ Franklin Court, Cambridge  [J Yang] 

    % initialise inputs 
    q_1_v = interpolate4q1 (wave,s_mech,options); 

    Phi_x = mk.Phi_x.d0;

    om_v = wave.F_f;

    aw = s_mech.waveH;

    g = options.g;
    wavetype = options.wavetype;
    Struct_type = options.Struct_type;

    N = s_geo_v.N;
    xi_v = s_geo_v.x;
    Li_v = s_geo_v.L;
    Di_v = s_geo_v.D(:,3);

    M = mk.M;
    K = mk.K;
    M_1 = mk.constrainM;
    K_1 = mk.constrainK;
    Ndof = mk.Ndof;

    Vb = mode.Vb;
    Nvb = mode.Nvb;
    D_om_b = mode.omb;

    Cm = s_mech.Morrison(2);
    Cd = s_mech.Morrison(3);
    rho_f = s_mech.rho(2);

    om_v = om_v(:)'; % make sure om_v is row vector

    Nom = numel(om_v);
    if Nom > 1
        dom = om_v(2)-om_v(1); % delta omega
    end

    Coef_iter = 0.5; % iteration steps <---------------------------
    Niter = 50;

    % define intial vectors
    k_w = om_v.^2/g;  % wavenumber (based on linear wave theory)

    u_i_v = om_v.*aw.*exp(-xi_v*k_w);  % fluid velocity, matrix format, N x Nom
    u_i_dot_v = 1i*om_v.^2.*aw.*exp(-xi_v*k_w); % fluid acceleration,  N x Nom

    u_i_v(xi_v<0,:) = 0;
    u_i_dot_v(xi_v<0,:) = 0;

    q_dot = zeros(Ndof,Nom); 

    if Struct_type == 1
        y_dot = Phi_x*q_dot;
    elseif Struct_type == 2
        y_dot = Phi_x(:,2:end)*q_dot+Phi_x(:,1)*(1i*om_v.*q_1_v);
    end
    ur = u_i_v - y_dot;  % relative velocity for ith strip

    % linearization coefficient
    if wavetype == 1
        gamma_old = 8/(3*pi)*abs(ur); % gamma is now in matrix form
    elseif wavetype == 2 
        % for random wave, gamma is over the whole frequency range
        Sxx = wave.Sxx.';
        Sur = abs(ur).^2.*Sxx;
        gamma_old = sqrt(8/pi)*sqrt(sum(Sur*dom,2)); % linearised coefficient, initialise (sqrt(8/pi)*std(ur))
    end

    % iteration to solve the equation because C is dependent on gamma
    for iter = 1:Niter

        H = zeros(N, Nom); %response spectrum,each row for each dof
        Hq = zeros(Ndof, Nom);

        ur = u_i_v - y_dot;  % relative velocity for ith strip

        if wavetype == 1

            gamma_new = 8/(3*pi)*abs(ur); % gamma is now in matrix form

        elseif wavetype == 2

            % for random wave, gamma is over the whole frequency range
            Sxx = wave.Sxx.';
            Sur = abs(ur).^2.*Sxx;
            gamma_new = sqrt(8/pi)*sqrt(sum(Sur*dom,2)); % linearised coefficient, initialise (sqrt(8/pi)*std(ur))

        end

        gamma_T = gamma_old + ...
            Coef_iter*(gamma_new-gamma_old); % iteration here

        if wavetype == 2 % for random wave, do this out of om loop to save time

            gamma = gamma_T;
            % Damping Matrix (needs update iteratively
            BDi_v = 1/2*rho_f*Di_v.*Li_v.*Cd.*gamma;    
            BDi_v(xi_v<0) = 0; % only effective below fluid surface 

            if Cd == 0  % if no drag force, add artificial damping
                C = 0.1*M;
            else
                % C=(Phi_x.')*diag(BDi_v)*Phi_x;
                C = (Phi_x.').*BDi_v.'*Phi_x;
            end

            if Struct_type == 2
                C_1 = C(2:end,1);
                C = C(2:end,2:end);
            end
        end

        % loop over omega
        jj = 0;

        for om = om_v % harmonic frequenccy from exp(i*om*t)

            jj = jj+1;

            q_1 = q_1_v(jj);

            if wavetype == 1
                gamma = gamma_T(:,jj); 
                % gamma_T is in matrix form for harmonic wave
                
                % Damping Matrix (needs update iteratively
                BDi_v = 1/2*rho_f*Di_v.*Li_v.*Cd.*gamma;
                BDi_v(xi_v<0) = 0; % only effective below fluid surface 


                if Cd == 0  % if no drag force, add artificial damping
                    C = 0.1*M;
                else
                    % C=(Phi_x.')*diag(BDi_v)*Phi_x;
                    C = (Phi_x.').*BDi_v.'*Phi_x;
                end

                if Struct_type == 2
                    C_1 = C(2:end,1);
                    C = C(2:end,2:end);
                end
            end

            u_v = u_i_v(:,jj);  % fluid velocity, vector of all strips
            u_dot_v = u_i_dot_v(:,jj); % fluid acceleration, vector of all strips

            % force vector (needs update iteratively
            FMi_v = pi/4*(Di_v).^2*rho_f*Cm.*Li_v.*u_dot_v +...
                1/2*rho_f*Cd.*Di_v.*Li_v.*u_v.*gamma;
            FMi_v(xi_v<0) = 0; % only effective below fluid surface

            F = (Phi_x.')*FMi_v;

            if Struct_type == 2 % for riser, there is presecried top motion

                F = F(2:end);

                q_1_2d = -om^2*q_1;
                q_1_1d = 1i*om*q_1;
                F = F - M_1*q_1_2d - K_1*q_1 - C_1*q_1_1d;%;

            end

            %----------------------------------------
            % calculate response at in modal coodinates
            Fb = Vb.'*F;  % force in modal coordinates
            Cb = Vb.'*C*Vb; % damping in modal coordinates

            D_b = -om^2*eye(Nvb) + 1i*om*Cb + diag(D_om_b);
            b_v = D_b\Fb;

            %----------------------------------------
            % mode acceleration
            q_v = K\(F-M*Vb*(-om^2*b_v)-C*Vb*(1i*om*b_v)); 
            % response in generalised coodindates
%             q_v=Vb*b_v; % if no mode acceleration


            q_dot_v = 1i*om*q_v; % generalised coor

            if Struct_type == 1

                y_dot(:,jj) = Phi_x*q_dot_v; % transform back to physical coor
                H(:,jj) = Phi_x*q_v; % y response (horizontal)

            elseif Struct_type == 2

                y_dot(:,jj) = Phi_x(:,2:end)*q_dot_v + Phi_x(:,1)*(1i*om*q_1);
                H(:,jj) = Phi_x(:,2:end)*q_v+Phi_x(:,1)*(q_1);
            end

            Hq(:,jj) = q_v; % y response (horizontal)
        end

        gamma_old = gamma_new; % save old gamma before next iteration cycle

    end

    if Struct_type == 2
        Hq = [q_1_v; Hq];
    end

 end