% Makes 36 plots of randomly sampled traces from the raw data before or
% after cleaning. Run after running CleanDVData.
% Paul Marchando
% 11-28-2018

which_plot = input('Do you want to view before or after cleaning? (b/a): ','s');

if which_plot == 'b'

indices_new = randperm(length(Fluo),36);

for i = 1:36
    
    subplot(6,6,i)
    plot(Fluo{indices_new(i)})
    hold on
    plot(smooth(Fluo{indices_new(i)}))
    fake_DV = num2str(DVpos{indices_new(i)});
    y_pos_str = num2str(ypos{indices_new(i)});
    title(['DV: ' fake_DV ' Y: ' y_pos_str])
    
end

else

indices_new = randperm(length(working_dat_final),36);

for i = 1:36
    
    subplot(6,6,i)
    plot(working_dat_final{indices_new(i),1})
    hold on
    plot(smooth(working_dat_final{indices_new(i),1}))
    fake_DV = num2str(working_dat_final{indices_new(i),2});
    y_pos_str = num2str(working_dat_final{indices_new(i),4});
    title(['DV: ' fake_DV ' Y: ' y_pos_str])
    
end

end