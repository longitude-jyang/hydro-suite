% callcode.m initialise the hydrodynamic calculation with simple options
% calCase:  three structures to chooose from
%        calCase ==1 --> a full scale rigid body floating wind turbine (tower only) 
%        calCase ==2 --> a model scale rigid body tensioned floating column (for wave tank testing)
%        calCase ==3 --> a full scale flexible tensioned marine riser 
% wavetype: choose either 
%        wavetype ==1 --> a harmonic regular wave 
%        wavetype ==2 --> a random wave (Jonswap spectrum) 
% calctype: choose either 
%        calcType ==1 --> frequency domain analysis with linearised equation 
%        calcType ==2 --> nonlinear time domain analysis by numerical integration

% 15/04/2021 @ Franklin Court, Cambridge  [J Yang] 

calCase  =  1  ;  % choose structure, FWT == 1, FWTtank == 2, riser == 3
wavetype =  1  ;  % choose wave type, regular == 1, random == 2
calcType =  1  ;  % choose freuqency == 1 or time == 2


 disp('                                                            ')
 disp(' --- Solving The Equations --- ')
 disp(' Analysis Starts: ...')
 tic;  
      %---------------------------------------------------------------------
      
       switch calCase 
            case 1
                ModPar = calFWT (wavetype, calcType); 
            case 2
                ModPar = calFWTtank (wavetype, calcType); 
            case 3
                ModPar = calRiser (wavetype, calcType); 
        end
        
        Out = maincode (ModPar);    
      %---------------------------------------------------------------------
     
 elapseTime = floor(toc*100)/100; Out.elapseTime =  elapseTime; 
 disp(strcat('Analysis Completed: ',num2str(elapseTime),'[s]'))
 
 % plot response, bending moment and bening stress 
 plotResponse; 