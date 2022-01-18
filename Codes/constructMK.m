function [mk, mode] = constructMK(s_geo, s_geo_v, s_mech, options)
% constructMK.m contruct the mass, stiffness matrix (MK)

% 07/07/2020 @ Franklin Court, Cambridge  [J Yang] 
% 09/04/2021 @ Franklin Court,, Cambridge  [J Yang]  --> code checked
% 18/01/2022 @ Franklin Court,, Cambridge  [J Yang]  --> use Bi_v, instead
% of Mia_v for Mi_t_3 

    % assign values 
    Struct_type = options.Struct_type;
    Ndof        = options.Ndof; % number of shape functions
    Nvb         = options.Nvb;  % number of modes to be kept
    g           = options.g; 


    L1          = s_geo.L1;
    L2          = s_geo.L2;
    L           = s_geo.L;

    N           = s_geo_v.N;
    xi_v        = s_geo_v.x;
    Li_v        = s_geo_v.L;
    
    xc          = s_geo_v.xc;
    xm          = s_geo_v.xm;

    Mi_v        = s_geo_v.m;
    Mia_v       = s_geo_v.ma;
    Bi_v        = s_geo_v.Bi;

    Ii_v        = s_geo_v.I;

    E           = s_mech.E;
    T0          = s_mech.T0;
    s_alpha     = s_mech.spring;

    %-------------------------------------------------------------------------
    % (1) Shape Functions   (case dependent)
    %-------------------------------------------------------------------------
   
    [Nphi,Phi_x] = constructShapeFunction (s_geo, s_geo_v, options);
    Phi_x_0d     = Phi_x.d0;
    Phi_x_1d     = Phi_x.d1;
    Phi_x_2d     = Phi_x.d2;

    %-------------------------------------------------------------------------
    % (2) Mass Matrix
    %-------------------------------------------------------------------------
    
    % summation in vector forms: NphixN * NxN * NxNphi--> NphixNphi
    M_T         = (Phi_x_0d.')*diag((Mi_v+Mia_v))*Phi_x_0d; 

    %-------------------------------------------------------------------------
    % (3) Stiffness Matrix
    %-------------------------------------------------------------------------

    % use this discretised version if non-unform
    Mi_t_3      = sum(triu(repmat((Mi_v-Bi_v),1,N))); 
    Mi_i2N      = diag(T0-Mi_t_3*g);

    %  bending contribution
    K1          = (Phi_x_2d.')*diag(E*Ii_v.*Li_v)*Phi_x_2d;  
    %  effective tension contribution
    K2          = (Phi_x_1d.')*(Mi_i2N)*(Phi_x_1d.*Li_v) ;  
    %  spring contribution
    K3          = s_alpha*(Phi_x_0d(xi_v==xm,:).')*Phi_x_0d(xi_v==xm,:) ; 

    K_T         = K1+K2+K3;

    
    %-------------------------------------------------------------------------
    % (4) Prescribed motion
    %-------------------------------------------------------------------------
    if Struct_type == 1  % FWT

        K = K_T;
        M = M_T;

        M_1 = [];
        K_1 = [];
        q_1 = 0;

    elseif Struct_type == 2  % riser

        if Nphi >= 2
            K = K_T(2:end,2:end); % contrained 1st rigid body mode
            M = M_T(2:end,2:end);

            M_1 = M_T(2:end,1);
            K_1 = K_T(2:end,1);        

            q_1 = 1;

        else
            M = M_T;
            K = K_T;

            M_1 = [];
            K_1 = [];
            q_1 = 0;
        end
    end

    %-------------------------------------------------------------------------
    % (5) Eigen Analysis
    %-------------------------------------------------------------------------
        [V,D_om] = eig(K,M);

        D_om = diag(D_om);
        [~,index] = sort(D_om);
        V = V(:,index);
        D_om = D_om(index);

        % Mass normalization
        for ii = 1 : Ndof
            V(:,ii) = V(:,ii)/sqrt((V(:,ii).'*M*V(:,ii)));
        end

    %-------------------------------------------------------------------------
    % (6) Truncate Modes
    %-------------------------------------------------------------------------
        Vb = V(:,1:Nvb); % truncated mode shapes

        % Mb=Vb.'*M*Vb;% should be unit matrix
        % Kb=Vb.'*K*Vb;  % should be D_om
        D_om_b = D_om(1:Nvb);

    %-------------------------------------------------------------------------
    % prepare for output
    %-------------------------------------------------------------------------
        mk.M = M;
        mk.K = K;
        mk.constrainM = M_1;
        mk.constrainK = K_1;
        mk.Phi_x = Phi_x;
        mk.Ndof = Ndof;
        mk.q_1 = q_1;

        mode.om = D_om;
        mode.omb = D_om_b;
        mode.Vb = Vb;
        mode.Nvb = Nvb;

end
    

