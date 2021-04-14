function BM = cal_BM_AC(s_geo_v, s_mech, mk, wave, response, options)
% cal_BM_AC calculates the bending moment along the structure
% two methods, one is EI based, when the 2nd derivative is avaiable
% 2ns method is IR based, where the moment is calculated along the
% structure by force balance, including interia froce (D'lembert)

% 15/03/2019 @ JDG-01, Cambridge  [J Yang] 
% 09/04/2021 @ Franklin Court,, Cambridge  [J Yang]  --> code checked


wavetype = options.wavetype;
TF = options.TF;
Nwin = options.Nwin;

E = s_mech.E;

Phi_x = mk.Phi_x;
Phi_x_1d = Phi_x.d1;
Phi_x_2d = Phi_x.d2;
Phi_x_3d = Phi_x.d3;
Phi_x_4d = Phi_x.d4;

%--------------------------------------------------------------------------
% (1)  get spatial derivatives of the response
%--------------------------------------------------------------------------
switch TF
    case 1
        Hq = response.Hq;
    case 2
        Hq = response.qf;
end

theta = Phi_x_1d*Hq;
theta2 = Phi_x_2d*Hq;
theta3 = Phi_x_3d*Hq;
theta4 = Phi_x_4d*Hq;

%--------------------------------------------------------------------------
% (2) calculates the bending moment/stress [M=EIy'']
%--------------------------------------------------------------------------

Ib = s_geo_v.I; % in vector form
Dout = s_geo_v.D(:,2);

BM_EI.Smm = -E*Ib.*theta2; % bending moment
BM_EI.Fs = E*Ib.*theta3; % shear force
BM_EI.Qs = E*Ib.*theta4; % external forces
BM_EI.Sigma = 1/2*E*Dout.*theta2*1e-6;  % bending stress  [N/mm2]

if TF==2 && wavetype==2
    
    om_range = wave.F_f;
    fs = om_range(end)/(2*pi)*2;
    
    Nt = length(response.tf);
    
    [Smm,~] = pwelch(BM_EI.Smm.',round(Nt/Nwin),[],[],fs);
    [Sigma,~] = pwelch(BM_EI.Sigma.',round(Nt/Nwin),[],[],fs);
    
    SmmEI = Smm.'/(2*pi);
    SigmaEI = Sigma.'/(2*pi);
    
end

%--------------------------------------------------------------------------
% (3) calculates the bending moment/stress [with Inertia Relief (IR)]
%--------------------------------------------------------------------------

BM_IR = cal_BM_IR(s_geo_v, s_mech, wave, theta, theta2, response, options);

BM_IR.Sigma = BM_IR.Smm.*1/2.*Dout./Ib*1e-6;

BM_IR.Sigma (isnan(BM_IR.Sigma)) = 0; % replace NaN with zero (for zero diameter dummy section, there might be NaN)

switch TF
    case 1
        switch wavetype
            
            case 1
                SmmIR = abs(BM_IR.Smm);
                SigmaIR = abs(BM_IR.Sigma);
                
                SmmEI = abs(BM_EI.Smm);
                SigmaEI = abs(BM_EI.Sigma);
                
            case 2
                
                Sxx = wave.Sxx;
                
                SmmIR = abs(BM_IR.Smm).^2.*Sxx.';
                SmmEI = abs(BM_EI.Smm).^2.*Sxx.';
                
                SigmaIR = abs(BM_IR.Sigma).^2.*Sxx.';
                SigmaEI = abs(BM_EI.Sigma).^2.*Sxx.';
                
                omega = wave.F_f;
                
                rms_SmmEI = sqrt(trapz(omega, SmmEI, 2));
                rms_SmmIR = sqrt(trapz(omega, SmmIR, 2));
                
                rms_SigmaEI = sqrt(trapz(omega, SigmaEI, 2));
                rms_SigmaIR = sqrt(trapz(omega, SigmaIR, 2));
        end
        
    case 2
        
        switch wavetype
            
            case 1
                
                SmmIR    = max(BM_IR.Smm,[],2); % get peak value for all columns
                SigmaIR  = max(BM_IR.Sigma,[],2);
                
                SmmEI    = max(BM_EI.Smm,[],2);
                SigmaEI  = max(BM_EI.Sigma,[],2);
                
            case 2
                
                [Smm,~]   = pwelch(BM_IR.Smm.',round(Nt/Nwin),[],[],fs) ;
                [Sigma,~] = pwelch(BM_IR.Sigma.',round(Nt/Nwin),[],[],fs) ;
                
                SmmIR     = Smm.'/(2*pi);
                SigmaIR   = Sigma.'/(2*pi);
                
                rms_SmmEI   = std(BM_EI.Smm.').';% std value for all strips
                rms_SmmIR   = std(BM_IR.Smm.').';
                rms_SigmaEI = std(BM_EI.Sigma.').';
                rms_SigmaIR = std(BM_IR.Sigma.').';
        end
end

%--------------------------------------------------------------------------
% (4) output results
%--------------------------------------------------------------------------

BM.SmmEI = SmmEI;
BM.SigmaEI = SigmaEI;

BM.SmmIR = SmmIR;
BM.SigmaIR = SigmaIR;

if wavetype==2
    
    BM.rms_SmmEI = rms_SmmEI;
    BM.rms_SmmIR = rms_SmmIR;
    BM.rms_SigmaEI = rms_SigmaEI;
    BM.rms_SigmaIR = rms_SigmaIR;
end

end