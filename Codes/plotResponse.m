% plotResponse.m displays the responses, including displacement, bending
% moment, bending stress, from maincode.m; 

% 06/04/2020 @ Franklin Court, Cambridge  [J Yang] 

BM       = Out.BM;
response = Out.response;
wave     = Out.wave;
mode     = Out.mode;
TF       = ModPar.calcType;
wavetype = ModPar.wavetype;

%--------------------------------------------------------------------------
% Display results  
%--------------------------------------------------------------------------
        x_coor.name  ='Frequency [rad/s]';
        x_coor.namef ='Frequency [Hz]';
        x_coor2.name ='Time [s]';

        y_coor.name  ='x [m]';
        
        
        
        z_coor.SmmEI = abs(BM.SmmEI); % bening moment
        z_coor.SmmIR = abs(BM.SmmIR);

        z_coor.SigmaEI = abs(BM.SigmaEI); % stress
        z_coor.SigmaIR = abs(BM.SigmaIR);
        
        z_coor.EpsilonIR = abs(BM.SigmaIR*1e6/ModPar.E); % strain 

        figprop.BMEIcbar ='EI [Nm]';   
        figprop.BMIRcbar ='IR [Nm]'; 
        
        figprop.SigmaEIcbar ='EI (Stress \sigma [N/mm2])';   
        figprop.SigmaIRcbar ='IR (Stress \sigma [N/mm2])'; 
        
        figprop.EpsilonEIcbar ='EI \epsilon [-]'; 
        figprop.EpsilonIRcbar ='IR \epsilon [-]'; 
        

        if TF == 1 

            x_coor.value = wave.F_f;
            y_coor.value = Out.xCoord;                   
           

            if wavetype == 1
                z_coor.value = abs(response.frf); 
               
                
                figprop.title ='Frequency Domain Analysis & Harmonic Wave';
                figprop.cbar ='y/a_w [m/m]';            
                
                figprop.BMtitle      ='FD - Harmonic - BM';
                figprop.Sigmatitle   ='FD - Harmonic - Stress';
                figprop.Epsilontitle ='FD - Harmonic - Strain';
            
            elseif wavetype == 2
                
                z_coor.value   = abs(response.Syy); 
                z_coor.rms_Syy = abs(response.rms_Syy); 
                
                z_coor.rms_SigmaEI = abs(BM.rms_SigmaEI);
                z_coor.rms_SigmaIR = abs(BM.rms_SigmaIR);

                
                figprop.title ='Frequency Domain Analysis & Random Wave';
                figprop.cbar  ='y [m]';  
                
                figprop.BMtitle      ='FD - Random - BM';
                figprop.Sigmatitle   ='FD - Random - Stress';
                figprop.Epsilontitle ='FD - Random - Strain';
      
            end


        elseif TF == 2 
            
            if wavetype == 1
                x_coor.value = wave.T_f;
                y_coor.value = Out.xCoord;
                z_coor.value = abs(response.y_peak); 

                figprop.title ='Time Domain Analysis & Harmonic Wave';
                figprop.cbar  ='y/a_w [m/m]';
                
                figprop.BMtitle      ='TD - Harmonic - BM';
                figprop.Sigmatitle   ='TD - Harmonic - Stress';
                figprop.Epsilontitle ='TD - Harmonic - Strain';
                
            elseif wavetype == 2
                
                x_coor.value   = response.fSyy.'*2*pi;
                y_coor.value   = Out.xCoord;
                z_coor.value   = response.Syy; 
                z_coor.rms_Syy = abs(response.rms_Syy); 

                
                z_coor.rms_SigmaEI = abs(BM.rms_SigmaEI);
                z_coor.rms_SigmaIR = abs(BM.rms_SigmaIR);

                figprop.title   ='Time Domain Analysis & Random Wave';
                figprop.cbar    ='y [m]';
                
                figprop.BMtitle      ='TD - Random - BM';
                figprop.Sigmatitle   ='TD - Random - Stress';
                figprop.Epsilontitle ='TD - Random - Strain';
            end

        end

        % call function to display the results 
        
        
        om_range = x_coor.value;
        [~,iom]  = min(abs(om_range-sqrt(mode.om(2))));
       
        figExtra.om  = om_range;
        figExtra.iom = iom;
        
        % displacement
        display_res(x_coor,y_coor,z_coor,figprop); 
        % bending moment and stress/strain
        display_BM(x_coor,y_coor,z_coor,figprop,figExtra,wavetype)