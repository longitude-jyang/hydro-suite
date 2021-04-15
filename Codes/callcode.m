 
 calCase  = 3; 

 disp('                                                            ')
 disp(' --- Solving The Equations --- ')
 disp(' Analysis Starts: ...')
 tic;  
      %---------------------------------------------------------------------
      
       switch calCase 
            case 1
                calFWT; 
            case 2
                calFWTtank; 
            case 3
                calRiser; 
        end
        
        Out = maincode (ModPar);    
      %---------------------------------------------------------------------
     
 elapseTime = floor(toc*100)/100; Out.elapseTime =  elapseTime; 
 disp(strcat('Analysis Completed: ',num2str(elapseTime),'[s]'))
 
 plotResponse; 