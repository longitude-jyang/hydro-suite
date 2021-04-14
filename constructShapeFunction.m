function [Nphi,Phi_x] = constructShapeFunction (s_geo, s_geo_v, options)
% constructShapeFunction prepares the shape function to transform
% coordinate from those along the structure to generalised coordinates 
% this shape functions are case dependent 

% 15/03/2019 @ JDG-01, Cambridge  [J Yang] 
% 09/04/2021 @ Franklin Court,, Cambridge  [J Yang]  --> update with shape
% function for FWTtank case 

Struct_type = options.Struct_type;
isModelTest = options.isModelTest;
bodyflex    = options.bodyflex;
Ndof        = options.Ndof;             % number of shape functions

L1          = s_geo.L1;
L2          = s_geo.L2;
L           = L1 + L2;

xi_v        = s_geo_v.x;
xc          = s_geo_v.xc;
N           = s_geo_v.N;

if Struct_type==1  % FWT
        Nphi=Ndof;

        if bodyflex==0 % rigid

           
            if isModelTest == 1 
                
                d           = s_geo.d;
                L_S         = s_geo.L(2);
                
                Lw = d - L_S; % wet section of the floating column 
                
                % shape function 
                Phi_x_0d_1 = [L_S*ones(nnz(xi_v < Lw),1); d - xi_v(xi_v >= Lw)]; 
                Phi_x_0d_1 (end) = 0; 
                Phi_x_0d_2 = [(d - L_S - xi_v(xi_v < Lw)) ; zeros(nnz(xi_v >= Lw),1)]; 
                
                Phi_x_0d=[Phi_x_0d_1  Phi_x_0d_2]; % size N x 2
                
                % first derivative of shape function 
                Phi_x_1d_1 = [zeros(nnz(xi_v < Lw),1) ; -ones(nnz(xi_v >= Lw),1)]; 
                Phi_x_1d_2 = [-ones(nnz(xi_v < Lw),1) ; zeros(nnz(xi_v >= Lw),1)] ; 
                
                Phi_x_1d=[Phi_x_1d_1 Phi_x_1d_2];
            
            else
                
                Phi_x_0d=[ones(N,1)/L (xi_v-xc)/L]; % size N x 2
                Phi_x_1d=[zeros(N,1) ones(N,1)/L];
                
            end
            
            Phi_x_2d=[zeros(N,1) zeros(N,1)];
            Phi_x_3d=[zeros(N,1) zeros(N,1)];
            Phi_x_4d=[zeros(N,1) zeros(N,1)];

        elseif bodyflex==1 % beam
            

            % cos(x) functions
            % shape function matrix, each column corresponds to one 'mode shape' for all coordinate of x
            lamda_v=[1:Nphi-2]*pi/L;

            Phi_x_0d=cos((xi_v+L1)*lamda_v); % form the shape function matrix: size N x Nphi
            Phi_x_0d=[ones(N,1)/L (xi_v+L1)/L  Phi_x_0d];

            Phi_x_1d=-repmat(lamda_v,N,1).*sin((xi_v+L1)*lamda_v); % 1st derivative of shape function
            Phi_x_1d=[zeros(N,1) ones(N,1)/L  Phi_x_1d];

            Phi_x_2d=repmat([ 0 0 -lamda_v.^2],N,1).*[ zeros(N,1) zeros(N,1) cos((xi_v+L1)*lamda_v)]; % 2nd derivative

            Phi_x_3d=-repmat([0 0 lamda_v.^3],N,1).*[ zeros(N,1) zeros(N,1) sin(xi_v*lamda_v)]; % 3rd derivative
            Phi_x_4d=-repmat([0 0 lamda_v.^4],N,1).*[zeros(N,1) zeros(N,1) cos(xi_v*lamda_v)]; % 4rd derivative

        end
       

    elseif Struct_type==2  % riser

        Nphi=Ndof+1;  % plus one, for rigid body motion

        % shape function matrix, each column corresponds to one 'mode shape' for all coordinate of x
        lamda_v=[1:Nphi-1]*pi/L;


        Phi_x_0d=sin(xi_v*lamda_v); % form the shape function matrix: size N x Nphi
%         Phi_x_0d(1,:) = 0;  
%         Phi_x_0d(end,:) = 0;        % due to the discrete coordinate, the shape function value at boundary are not zero
        Phi_x_0d=[1-xi_v/L  Phi_x_0d];

        Phi_x_1d=repmat(lamda_v,N,1).*cos(xi_v*lamda_v); % 1st derivative of shape function
        Phi_x_1d=[-ones(N,1)/L  Phi_x_1d];

        Phi_x_2d=-repmat([0 lamda_v.^2],N,1).*[ zeros(N,1) sin(xi_v*lamda_v)]; % 2nd derivative
        Phi_x_3d=-repmat([0 lamda_v.^3],N,1).*[ zeros(N,1) cos(xi_v*lamda_v)]; % 3rd derivative
        Phi_x_4d=repmat([0 lamda_v.^4],N,1).*[zeros(N,1) sin(xi_v*lamda_v)]; % 3rd derivative



end
    
Phi_x=struct('d0',{Phi_x_0d},...
    'd1',{Phi_x_1d},...
    'd2',{Phi_x_2d},...
    'd3',{Phi_x_3d},...
    'd4',{Phi_x_4d});