function plotMultipleTabsDorsal(dataType1, dataType2)

addDVStuffToSchnitzCells(dataType1)
addDVStuffToSchnitzCells(dataType2)

[~,~,axs0, axs0mrna] = plotFracByDlFluo(dataType1);
[~,~, axs1, axs1mrna] = plotFracByDlFluo(dataType2);

copyPlot(axs0{1}, axs1{1});
copyPlot(axs0{2}, axs1{2});
legend(axs1{1}, '1DG', '0DG');


copyPlot(axs0mrna{1}, axs1mrna{1});
copyPlot(axs0mrna{2}, axs1mrna{2});
legend(axs1mrna{1}, '1DG', '0DG');

end
