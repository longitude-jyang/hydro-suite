function [x, t] = psd2tim(Sxx, om_v, Phase)
% reconstruct one time series data from given PSD spectrum 
% INPUTS
%       Sxx - One-sided energy spectral density [meter^2*sec]
%       om_v -  frequency spectrum [rad/sec]
%       Phase - random phase [rad]
%
% OUTPUTS
%       x - Generated timeseries
%       t - Time vector corresponding to timeseries generated

% 21/03/2019 @ JDG-01, Cambridge  [J Yang] 


    Sxx2=[Sxx ; repmat(zeros(size(Sxx)),1,1)]; % padding with zeros
    Phase2=[Phase ;repmat(zeros(size(Sxx)),1,1)];
    
    [N,~]= size(Sxx2);    
    dom = om_v(2)-om_v(1);           % Delta frequency [rad/sec]
    dt = 2*pi/(N*dom);               % Delta time [sec]

    t = (0:N-1)*dt;                  % Time vector
    T = (t(N)-t(1)).*N./(N-1);       % length of the output timeseries [seconds]
    X_om =sqrt(2*Sxx2*dom);           % Magnitude of Linear Spectrum [meters]

    
    t=t.'; % colum vector 
    x = N*real(ifft(X_om.*exp(1i*Phase2))); % time series, with ifft (ifft normalied by N in Matlab, so df*N=1/dt) 

%     % ----------------------------------------------------------------------
%     % lines below are commennted off, but can be switched back on to check by converting back to PSD 
%     %%
%         A=1/N*fft(x);
%         Axx=1/(dom)* abs(A).^2;
%     
%         df=1/T;
%         f=[0:df:(N-1)*df]';
%     
%         N2=ceil(N/2);
%         f2=f(1:N2);
%         Axx=2*Axx(1:N2,:);
%    
%         figure
%         subplot(211)
%         plot(om_v, Sxx(:,end),2*pi*f2,Axx(:,end),'r--')
%         legend('PSD','Recovered-PSD')
%         xlabel('Frequency (rad/s)')
%         ylabel('Sxx [m^2s]')
% 
%         subplot(212)
%         plot(t, x)
%         title('Generated Timeseries From PSD')
%         xlabel('Time (Sec)')
%         ylabel('Amplitude [m]')