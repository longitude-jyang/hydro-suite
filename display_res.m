% diplay the results

function display_res(x_coor,y_coor,z_coor,figprop)

x=x_coor.value;
y=y_coor.value;
z=z_coor.value;

figure 
     
surf(x,y,z,'EdgeColor','none');   
axis xy; axis tight; colormap(jet); view(0,90);
set(gca, 'YDir','reverse')
title (figprop.title)

xlabel(x_coor.name)
ylabel(y_coor.name)


shading interp
hcb=colorbar;
title(hcb,figprop.cbar)
set(gca,'FontSize',14)