% diplay the results for bending moment and bending stress/strain

function display_BM(x_coor,y_coor,z_coor,figprop,figExtra,wavetype)


    x = x_coor.value;
    y = y_coor.value;

    SmmEI   = z_coor.SmmEI;
    SmmIR   = z_coor.SmmIR;
    SigmaEI = z_coor.SigmaEI;
    SigmaIR = z_coor.SigmaIR;

    epsilonIR = z_coor.EpsilonIR;

    if wavetype == 1
        dis = z_coor.value;
    elseif wavetype == 2
        dis = z_coor.rms_Syy;
    end

%-----------------------------
% 2-D plot 

    % compare bending moment (new figure)
    compare2Dsurf (x_coor,y_coor,SmmEI,SmmIR,figprop.BMtitle,figprop.BMtitle,figprop.BMEIcbar,figprop.BMIRcbar) ; 
    
    % compare bending stress (new figure)
    compare2Dsurf (x_coor,y_coor,SigmaEI,SigmaIR,figprop.Sigmatitle,figprop.Sigmatitle,figprop.SigmaEIcbar,figprop.SigmaIRcbar) ; 

    
%-----------------------------
% 1-D plot at specific frequency (one of the natural frequency)
    
    iom=figExtra.iom;
    om_range=figExtra.om;
    
    % for regular wave
    if wavetype==1  


        figure
        subplot(121)
        plot(y,dis(:,iom),'b')

        xlabel(y_coor.name)
        ylabel('Displacement [m]')
        title([{figprop.title} ; ...
            {strcat('Displacement for Frequency @',num2str(om_range(iom)))}])

        set(gca, 'XDir','reverse')
        set(gca,'view',[90 -90])
        set(gca,'FontSize',14)

        subplot(122)
        plot(y,SigmaEI(:,iom),'r--',...
            y,SigmaIR(:,iom),'k')

        legend('EI','IR')

        xlabel(y_coor.name)
        ylabel('Stress \sigma [N/mm2]')
        title(strcat('\sigma for Frequency @',num2str(om_range(iom))))

        set(gca, 'XDir','reverse')
        set(gca,'view',[90 -90])
        set(gca,'FontSize',14)

    elseif wavetype==2

        figure
        plot(y,SigmaEI(:,iom),'r--',...
            y,SigmaIR(:,iom),'k')

        legend('EI','IR')

        xlabel(y_coor.name)
        ylabel('Stress \sigma [N/mm2]')
        title(strcat('\sigma for Frequency @',num2str(om_range(iom))))

        set(gca, 'XDir','reverse')
        set(gca,'view',[90 -90])
        set(gca,'FontSize',14)


        figure
        subplot(121)
        plot(y,dis,'b')

        xlabel(y_coor.name)
        ylabel('Displacement [m]')
        title([{figprop.title} ; ...
            {'r.m.s displacement'}])

        set(gca, 'XDir','reverse')
        set(gca,'view',[90 -90])
        set(gca,'FontSize',14)

        subplot(122)
        plot(y,z_coor.rms_SigmaEI,'r--',...
            y,z_coor.rms_SigmaIR,'k')
        legend('EI','IR')


        xlabel(y_coor.name)
        ylabel('Stress [N/mm2]')
        title('r.m.s. Stress')

        set(gca, 'XDir','reverse')
        set(gca,'view',[90 -90])
        set(gca,'FontSize',14)

    end
    
end

function compare2Dsurf (x_coor,y_coor,S1,S2,figtitle1,figtitle2,bartitle1,bartitle2)    
    x = x_coor.value;
    y = y_coor.value;


    figure 

    subplot(211)
    surf(x,y,S1,'EdgeColor','none');   
        axis xy; axis tight; colormap(jet); view(0,90);
        set(gca, 'YDir','reverse')
        title (figtitle1)

        xlabel(x_coor.name)
        ylabel(y_coor.name)


        shading interp
        hcb=colorbar;
        title(hcb,bartitle1)
        set(gca,'FontSize',14)

    subplot(212)
    surf(x,y,S2,'EdgeColor','none');   
        axis xy; axis tight; colormap(jet); view(0,90);
        set(gca, 'YDir','reverse')
        title (figtitle2)

        xlabel(x_coor.name)
        ylabel(y_coor.name)


        shading interp
        hcb=colorbar;
        title(hcb,bartitle2)
        set(gca,'FontSize',14)
end



