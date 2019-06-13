% Makes 16 plots of randomly sampled Gaussian fits. Run after running
% DVpolyfit.
% Paul Marchando
% 11-28-2018

ind = randperm(length(gauss_out),16);

for q = 1:16
    
    temp_range = find(fake_dv_comb(:,3) == ind(q));
    subplot(4,4,q)
    scatter(fake_dv_comb(temp_range(1):temp_range(end),1),fake_dv_comb(temp_range(1):temp_range(end),2))
    hold on
    plot(gauss_out{ind(q),1});
    title(['Frame: ' num2str(ind(q))])
    
end