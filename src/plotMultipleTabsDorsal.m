function plotMultipleTabsDorsal(dataType1, dataType2)
% 
addDVStuffToSchnitzCells(dataType1)
addDVStuffToSchnitzCells(dataType2)

[~,~,axs1, axs1mrna] = plotFracByDlFluo(dataType1);
[~,~, axs2, axs2mrna] = plotFracByDlFluo(dataType2);

copyPlot(axs1{1}, axs2{1});
copyPlot(axs1{2}, axs2{2});
legend(axs2{1}, dataType2, dataType1);


copyPlot(axs1mrna{1}, axs2mrna{1});
copyPlot(axs1mrna{2}, axs2mrna{2});
legend(axs2mrna{1}, dataType2, dataType1);

end
