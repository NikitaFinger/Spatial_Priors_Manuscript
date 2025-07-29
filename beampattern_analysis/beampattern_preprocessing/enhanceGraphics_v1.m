function enhanceGraphics(outputDir, out_fname, frames, bat, bpp, side_view, aa)
    % Create a new figure for enhanced graphics
    figure; 
    set(gcf, 'pos', [10 40 1200 800], 'color', 'w');
    set(gca, 'position', [.07 .06 .88 .89], 'units', 'normalized');
    hold on;

    % Plot microphones
    goodch = 1:length(bpp.mic_loc(:,1)); % Alternatively load in the good ch from the ch_ex var
    plot3(bpp.mic_loc(goodch, 1), bpp.mic_loc(goodch, 2), bpp.mic_loc(goodch, 3), ...
        'o', 'MarkerSize', 8, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'r');

    % Plot bat's path
    plot3(bat(frames(1):frames(end), 1), bat(frames(1):frames(end), 2), bat(frames(1):frames(end), 3), ...
        '-', 'linewidth', 3, 'color', 'k');

    if side_view
        view(0, 0)
    else
        view(2)
    end

    axis equal, grid on;
    set(gca, 'fontsize', 20, ...
        'XTick', aa(1):.4:aa(2), 'xticklabel', (aa(1):.4:aa(2))-aa(1), ...
        'YTick', aa(3):.4:aa(4), 'yticklabel', (aa(3):.4:aa(4))-aa(3), ...
        'GridLineStyle', '--', 'GridColor', 'k', 'GridAlpha', 0.5);
    title('Bat Flight Path and Beam Pattern', 'fontsize', 24, 'FontWeight', 'bold');
    xlabel('X Position (m)', 'FontSize', 20);
    ylabel('Y Position (m)', 'FontSize', 20);
    zlabel('Z Position (m)', 'FontSize', 20);
    box on;

    % Save the enhanced figure as an image or a new video frame
    saveas(gcf, fullfile(outputDir, [out_fname '_enhanced.png']));
    close(gcf);
end
