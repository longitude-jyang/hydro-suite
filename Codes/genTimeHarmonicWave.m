function [u_v,u_v_dot] = genTimeHarmonicWave (aw,om,t,phi,xi_v,g) 

    k_w     = om^2/g;  % wavenumber (based on linear wave theory)

    % fluid wave velocity/acceleration, vector of all strips
    u_v     = om*aw*exp(-k_w*xi_v)*cos(om*t + phi);  % fluid wave velocity, vector of all strips
    u_v_dot = -om^2*aw*exp(-k_w*xi_v)*sin(om*t + phi);  % wave acceleration 

    % above water, switch off fluid data
    u_v(xi_v<0,:)     = 0;
    u_v_dot(xi_v<0,:) = 0;
    
    
