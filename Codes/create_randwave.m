function [u, u_dot, t, Sxx, Suu, Sudud, om_v, TLP_R] = create_randwave(xi_v, s_mech, options)
% create_randwave create random wave in both frequency and time domain 
% at the moment, only Jonswap wave spectrum implemented, but additional
% types of sea states can be easily added. 

% 15/03/2019 @ JDG-01, Cambridge  [J Yang] 
% 09/04/2021 @ Franklin Court,, Cambridge  [J Yang]  --> code checked


    wavecode = options.wavecode;
    g = options.g;
    Nfft = options.Nfft;

    JONSWAP_w = s_mech.waveR; % for random wave

    xi_v = xi_v(:)'; % make sure xi_v is row vector

    if wavecode==1

        Hs = JONSWAP_w.Hs;

        % Two jonswap cases are considered in the code. If Tz is defined, the
        % jonswap spectrum used is the one defined in the DNVGL report. 

        if isfield(JONSWAP_w, 'Tz')

            Tz = JONSWAP_w.Tz;

            % Relation from DNVGL report  
            Tp = 1.2859*Tz;        
            wp = 2*pi/Tp;

            % Sampling performed up to 5 times the peak frequency
            wMax = options.wLim*wp;
            om_v = linspace(0, wMax, Nfft).';

            % Subfunction with the two parameters jonswap expression is called
            [Sxx, ~, Phase] = jonswapWithHs_AC(om_v, wp, Hs,options);

        else

            % This jonswap expression assumes a relation between Hs and wp 
            wp = 1.416/sqrt(Hs);

            % Sampling performed up to 5 times the peak frequency
            wMax = 5*wp;
            om_v = linspace(0, wMax, Nfft).';

            % Subfunction with the one parameters jonswap expression is called
            [Sxx, ~, Phase] = jonswap_AC(om_v, wp,options);

        end
    end

    % prescribed motion needs to be interpolated
    riser = s_mech.riser;
    H_TLP_v = riser.H_TLP; % frequency dependent
    T_TLP_v = riser.T_TLP;
    gamma_TLP_v = riser.gamma_TLP;

    H_TLP = interp1(2*pi./T_TLP_v, H_TLP_v, om_v, 'linear', 'extrap'); % column vector
    gamma_TLP = interp1(2*pi./T_TLP_v, gamma_TLP_v, om_v, 'linear', 'extrap'); % column vector

    
    
    % create wave spectrum
    k_w = om_v.^2/g; % wavenumber

    % do vector .* matrix which is vector repeated * matrix (each column in this case)
    % the resulted matrix has each colum represent a spectrum for each xi
    Suu = om_v.^2.*Sxx.*exp(-2*k_w*xi_v); % wave velocity spectrum
    Sudud = om_v.^4.*Sxx.*exp(-2*k_w*xi_v); % wave acceleration spectrum

        % no wave above the fluid surface
        Suu(:,xi_v<0)   = 0;
        Sudud(:,xi_v<0) = 0;
        
        
    % convert spectrum to random signals in the time domain 
    [u, t]     = psd2tim(Suu, om_v, repmat(Phase,1,length(xi_v)));
    [u_dot, ~] = psd2tim(Sudud, om_v, repmat(Phase+pi/2,1,length(xi_v)));

    [q_1, ~]    = psd2tim(Sxx.*H_TLP.^2, om_v, Phase-gamma_TLP);
    [q_1_1d, ~] = psd2tim(om_v.^2.*Sxx.*H_TLP.^2, om_v, Phase+pi/2-gamma_TLP);
    [q_1_2d, ~] = psd2tim(om_v.^4.*Sxx.*H_TLP.^2, om_v, Phase+pi-gamma_TLP);

    TLP_R.q_1 = q_1;
    TLP_R.q_1_1d = q_1_1d;
    TLP_R.q_1_2d = q_1_2d;

end

% jonswap_AC code is based on jonswap code but does not use the exact same
% expression 
function [S, Amp, Phase] = jonswap_AC(om_v, wm,options)

    g = options.g;
    alpha = 0.0081;
    gamma = 3.3;

    N = length(om_v);
    S = zeros(N, 1);

    for iW = 1:length(om_v)

        thisW = om_v(iW);

        if thisW <= wm
            sig = 0.07;
        else
            sig = 0.09;
        end

        num = (thisW - wm)^2;
        den = 2*sig^2*wm^2;

        gammaExp = gamma^exp(-num/den);

        if thisW==0
            S(iW)=0;
        else

            factorExp = (5/4)*(wm/thisW)^4;

            S(iW) = (alpha*g^2/thisW^5)*exp(-factorExp)*gammaExp;

        end
    end

    % Determine the frequency step from the frequency vector. Note that the
    % highest frequency step is extrapolated.
    domg = zeros(N,1);
    domg(1:end-1) = diff(om_v);
    domg(end) = domg(end-1);

    % Determine the amplitudes from the spectral values
    Amp = sqrt(2*S.*domg);

    % Random phases
    Phase = rand(length(om_v),1)*2*pi;

end

% jonswapWithHs_AC code is based on jonswap code but uses a two parameters
% expression for the spectrum
function [S, Amp, Phase] = jonswapWithHs_AC(om_v, wp, Hs,options)
    g = options.g;

    N = length(om_v);
    S = zeros(N, 1);

    for iW = 1:length(om_v)

        thisW = om_v(iW);

        gamma = 3.3;

        if thisW <= wp
            sig = 0.07;
        else
            sig = 0.09;
        end    

        if thisW==0
            S(iW)=0;
        else

            aGamma = (1 - 0.287*log(gamma));

            aPM = (5/16)*Hs^2*wp^4;

            factorExp1 = (5/4)*(wp/thisW)^4;

            factorExp2 = -0.5*(thisW - wp)^2/(sig*wp)^2;

            S(iW) = aGamma*aPM*thisW^(-5)*exp(-factorExp1)*gamma^exp(factorExp2);

        end
    end

    % Determine the frequency step from the frequency vector. Note that the
    % highest frequency step is extrapolated.
    domg = zeros(N,1);
    domg(1:end-1) = diff(om_v);
    domg(end) = domg(end-1);

    % Determine the amplitudes from the spectral values
    Amp = sqrt(2*S.*domg);

    % Random phases
    Phase = rand(length(om_v),1)*2*pi;

end



