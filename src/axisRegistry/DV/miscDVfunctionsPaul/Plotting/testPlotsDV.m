% Makes 36 plots of randomly sampled traces from the cleaned data with
% polynomial fits applied and fluorescence values determined. Run after
% running DVpolyfit.
% Paul Marchando
% 11-28-2018

indices = randperm(length(working_dat_final),36);

for i = 1:36
    
    subplot(6,6,i)
    plot(polyfit_data(indices(i),:))
    hold on
    plot(polyval(polyfit_out(indices(i),:),1:length(polyfit_data(1,:))))
    hold on
    plot(polyfit_maximums(indices(i),1),polyfit_maximums(indices(i),2),'r*');
    fake_DV = num2str(working_dat_final{indices(i),2});
    y_pos_str = num2str(working_dat_final{indices(i),4});
    title(['DV: ' fake_DV ' Y: ' y_pos_str])
    
end
