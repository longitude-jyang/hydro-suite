function  wave = constructWave (s_geo_v, s_mech, options)
% constructWave.m contruct ocean wave environment 

% 07/07/2020 @ Franklin Court, Cambridge  [J Yang] 

    wavetype = options.wavetype;

    switch wavetype

        % regular/harmonic wave
        case 1
            wave = constructWaveRegular (s_mech, options) ;

        % random wave
        case 2

            wave = constructWaveRandom (s_geo_v, s_mech, options) ;

    end
end

function wave = constructWaveRegular (s_mech, options) 

        om_range         = options.om_range;
        tmax             = options.tmax; 
        dtScale           = options.dtScale; 

        % for regular wave, frequency range needs to be defined for the
        % analysis (the interest could be around the fundamental mode of the structure)
        om_range_f       = om_range;
        N_om_f           = numel(om_range_f);
        
        % have differnet resolution for time analysis because time analysis would take very long time
        if numel(om_range) > N_om_f 
            om_range_t       = linspace(min(om_range),max(om_range),N_om_f);
        else
            om_range_t       = om_range; 
        end
        N_om_f = numel(om_range_t);
        
        % tmax should be at least 4 times of period of lowest frequency
        tspan            = 0:2*pi/max(om_range_t)/dtScale:tmax;
        Nt = numel(tspan);
        
        % prescribed motion needs to be interpolated
        aw               = s_mech.waveH;
        riser            = s_mech.riser;
        
        % frequency dependent TLP
        H_TLP_v          = riser.H_TLP; 
        T_TLP_v          = riser.T_TLP;
        gamma_TLP_v      = riser.gamma_TLP;
        
        % column vector
        H_TLP            = interp1(2*pi./T_TLP_v,H_TLP_v, ... 
                                   om_range_t,'linear','extrap'); 
        gamma_TLP        = interp1(2*pi./T_TLP_v,gamma_TLP_v, ... 
                                   om_range_t,'linear','extrap');
        
        % prescribed motion, size of Nom x Nt
        q_1              = aw*H_TLP.'.* ...
                           cos(om_range_t.'*tspan - ...
                           repmat(gamma_TLP.',1,Nt)); % Nom x Nt
        q_1_1d           = - om_range_t.'.*aw.*H_TLP.'.* ... 
                           sin(om_range_t.'*tspan - ...
                           repmat(gamma_TLP.',1,Nt));
        q_1_2d           = - om_range_t.'.^2.*aw.*H_TLP.'.* ...
                           cos(om_range_t.'*tspan - ...
                           repmat(gamma_TLP.',1,Nt));
        
        %--------------------------------------------------------------------------
        % prepare for output
        %--------------------------------------------------------------------------
        wave.Name        = 'Harmonic';
        
        wave.F_f         = om_range_f;
        wave.Nf          = N_om_f;
        
        wave.T_f         = om_range_t;
        wave.t           = tspan;
        wave.tmax        = tmax;
        wave.Nt          = Nt;
        
        TLP_H.q_1        = q_1;
        TLP_H.q_1_1d     = q_1_1d;
        TLP_H.q_1_2d     = q_1_2d;
        wave.TLP         = TLP_H;
end

function wave = constructWaveRandom (s_geo_v, s_mech, options) 
     
        % choose standard wave spectrum for random wave (the frequency range 
        % is based on the wavecode) JONSWAP wave
        xi_v             = s_geo_v.x;
        g                = options.g;
        
        % create random waves
        % use sub-routine to create wave time series
        % u & u_dot has each colum represents the time series of each strip
        % Sxx is the random wave spectrum
        [u_v, u_v_dot, tw, Sxx, Suu, Sudud, om_range, TLP_R] = ...
            create_randwave(xi_v, s_mech, options);  

        om_range        = om_range.'; % turn to row vector
        Nt              = length(tw);        
        k_w             = om_range.^2/g;
        
        wave.Name       = 'Random';        
        wave.F_f        = om_range;
        wave.k          = k_w;
        
        wave.t          = tw;
        wave.Nt         = Nt;
        
        wave.u          = u_v;
        wave.u_d        = u_v_dot;
        
        wave.Sxx        = Sxx;
        wave.Suu        = Suu;
        wave.Sudud      = Sudud;
        
        wave.TLP        = TLP_R;
end