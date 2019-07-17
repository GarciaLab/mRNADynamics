function plotMultipleTabsDorsal(dataType1, dataType2)

[~,~,axs0] = plotFracByDlFluo(dataType1);
[~,~, axs1] = plotFracByDlFluo(dataType2);

copyPlot(axs0{1}, axs1{1});
copyPlot(axs0{2}, axs1{2});
legend(axs1{1}, '1DG', '0DG');

end
