function Texture_plot(Nx,Ny,Sx,Sy,Sz)
    %This function plots a given spin texture
    %as a quiver plot of the x and y components and
    %a contour plot of the z component 

    %Defines the lattice
    [XX,YY] = meshgrid(0:1:Nx-1,0:1:Ny-1);
    
    %Plotting parameters
    scrsz = get(0,'ScreenSize');    
    figure('Position',[50 scrsz(4)/2+25 0.8*scrsz(3)/2 0.8*scrsz(4)/2])
    set(gca,'nextplot','replacechildren');
    set(gcf,'Renderer','zbuffer');
    
    %pcolor(X,Y,Sz); colorbar
    %shading('interp')
    %colormap(summer);

    %Plotting
    contourf(XX,YY,Sz,30,'Linestyle','None'); 
    c = colorbar;
    colormap(summer(80));
    clim([-1 1]);
    hold on

    %quiver(XX,YY,Sx,Sy,0,'-k');
    quiver(XX,YY,Sx,Sy,0.5,'-k')
    
    %Size parameters
    axis equal
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    c.FontSize = 12;
    hold off


end