

figs = plotDorsalDynamics('1Dg');
figs = plotDorsalDynamics('1Dg_2xDl', 'figs', figs);

for i = 1:length(figs)
    figure(figs{i})
    [~,icons,~,~] =  legend('1x', '2x');
    icons(3).Color = 'b';
    icons(5).Color = 'r';
end
