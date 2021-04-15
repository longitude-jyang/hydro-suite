% interpolate4q1.m interpolates the q1 motion for riser, with prescribed
% motion of top platform 

% 23/07/2020 @ Franklin Court, Cambridge  [J Yang] 

function q_1 = interpolate4q1 (wave,s_mech,options)
        om_range = wave.F_f;
        aw = s_mech.waveH;

        
        % prescribed motion needs to be interpolated
        riser = s_mech.riser;
        H_TLP_v = riser.H_TLP; % frequency dependent
        T_TLP_v = riser.T_TLP;
        gamma_TLP_v = riser.gamma_TLP;
        
        switch options.typeIntTLP
            
            case 'no extrapolation'
                H_TLP = interp1(2*pi./T_TLP_v, H_TLP_v, om_range, 'linear'); % column vector
                gamma_TLP = interp1(2*pi./T_TLP_v, gamma_TLP_v, om_range, 'linear'); % column vector
                
                H_TLP(isnan(H_TLP)) = 0;
                gamma_TLP(isnan(gamma_TLP)) = 0;
                
            case 'decay'

                H_TLP_v = H_TLP_v(2:end);
                T_TLP_v = T_TLP_v(2:end);
                gamma_TLP_v = gamma_TLP_v(2:end);
                
                H_TLP = interp1(2*pi./T_TLP_v, H_TLP_v, om_range, 'linear', 'extrap'); % column vector
                gamma_TLP = interp1(2*pi./T_TLP_v, gamma_TLP_v, om_range, 'linear', 'extrap'); % column vector

                H_TLP(H_TLP<0) = 0;

            case 'nearest'
                
                H_TLP = interp1(2*pi./T_TLP_v, H_TLP_v, om_range, 'nearest', 'extrap'); % column vector
                gamma_TLP = interp1(2*pi./T_TLP_v, gamma_TLP_v, om_range, 'nearest', 'extrap'); % column vector
                
            case 'decay with plateau'
                
                H_TLP_v = H_TLP_v(2:end);
                T_TLP_v = T_TLP_v(2:end);
                gamma_TLP_v = gamma_TLP_v(2:end);
                
                H_TLP = interp1(2*pi./T_TLP_v, H_TLP_v, om_range, 'linear', 'extrap'); % column vector
                gamma_TLP = interp1(2*pi./T_TLP_v, gamma_TLP_v, om_range, 'linear', 'extrap'); % column vector

                H_TLP(H_TLP<0) = 0;
                H_TLP(H_TLP> max(H_TLP_v)) = max(H_TLP_v);
                    
        end
        
        q_1 = aw*H_TLP.*exp(-1i*gamma_TLP); % frequency dependent