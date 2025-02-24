function plot_nelder_mead_3d(f, xlim, ylim, zlim)

    % xlim and ylim are the limits of the plot for x and y axis
    [X, Y] = meshgrid(linspace(xlim(1), xlim(2), 100), linspace(ylim(1), ylim(2), 100));
    Z = arrayfun(@(x, y) f([x; y]), X, Y);

    figure;
    hold on;
    surf(X, Y, Z, 'FaceAlpha', 0.7, 'EdgeColor', 'none');
    colormap turbo;
    colorbar;
    title('3D Visualization of Nelder-Mead Optimization');
    xlabel('x_1');
    ylabel('x_2');
    zlabel('f(x)');
    view(135, 30);
    shading interp;
end
