dataTypes = {'1Dg_2xDl', '1DgW_and_1DgW2xDl', '1DgVW_FFF'}; %affinities
% dataTypes = {'1Dg', '1Dg-5_FFF'};
% dataTypes = {'1Dg-5_FFF', '1Dg_2xDl', '1DgW_and_1DgW2xDl', '1DgVW_FFF'};
activities = {'fraction', 'max', 'timeon'};


clrmp = single(lines(length(dataTypes)));

for j = 1:length(activities)
    for i = 1:length(dataTypes)
        if j == 1 & i > 2
%             compileAllProjects(dataTypes{i})
% %             binDorsal(dataTypes{i}, true)
        end
        plotFracByDlFluo2(dataTypes{i}, activities{j});
        if i == 1
            ax1 =gca;
            for k = 1:length(ax1.Children)
%                 set(ax1.Children, 'Color', clrmp(1), 'MarkerFaceColor', clrmp(1));
                if isprop(ax1.Children(k), 'Color')
                    set(ax1.Children(k), 'Color', clrmp(1, :));
                end
                if isprop(ax1.Children(k), 'MarkerFaceColor')
                    set(ax1.Children(k), 'MarkerFaceColor', clrmp(1, :));
                end
            end
        else
            ax = gca;
            copyPlot(ax, ax1, 'close', 'colorMap', clrmp(i, :));
            xlim([0, 4000]);
            ymax = get(gca, 'YLim');
            ylim([0, ymax(2)]);
            leg = get(gca, 'Legend'); w=.02;h=.01;set(leg, 'Units', 'normalized', 'Position', [1-w, 1-h,w, h], 'Box','off');
            set(leg, 'Visible', 'off');
        end
    end
end
