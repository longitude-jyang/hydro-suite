function response = solve4Response(s_geo, s_geo_v, s_mech, wave, mk, mode, options)
% solve4Response.m solve the structural response due to wave exctations, in
% either time or frequency domain

% 23/07/2020 @ Franklin Court, Cambridge  [J Yang] 

    wavetype = options.wavetype;
    TF       = options.TF;
    

    switch TF
        case 1  % choose Frequency domain analysis
            
            % for frequency analysis, solve for FRF first, regardless of
            % wave types 
            [Hfrf, Hq] = solve4frf(s_geo_v, s_mech, mk, mode, wave, options);
            
            response.frf  = Hfrf;
            response.Hq   = Hq;
            
            switch wavetype

                % regular/harmonic wave
                case 1
                   


                % random wave
                case 2

                    om_range = wave.F_f;
                    
                    Sxx = wave.Sxx;                % wave spectrum 
                    Syy = abs(Hfrf).^2.*Sxx.';     % response spectrum
                    response.Syy = Syy;

                    rms_Syy = sqrt(sum(Syy*(om_range(2)-om_range(1)),2)); % std value for all strips
                    response.rms_Syy = rms_Syy;
                                

            end
            
        case 2 % choose Time domain analysis
            
            switch wavetype

                % regular/harmonic wave
                case 1
                    
                    % for harmonic wave, the time response is computed for
                    % each frequency interested in time domain, so the
                    % bending moment/stress is also computed here 
                    response = solve4HarmonicTimeResponse( s_geo_v, s_mech, wave, mk, mode, options);

                % random wave
                case 2
                    
                    response = solve4RandomTimeResponse(s_geo_v, s_mech, wave, mk, mode, options);


            end
            
            
    end
end