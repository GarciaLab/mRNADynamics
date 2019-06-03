for q = 1:length(gauss_out)
    
    temp_range = find(fake_dv_comb(:,3) == q);
    if isempty(temp_range) == 0
    scatter(fake_dv_comb(temp_range(1):temp_range(end),1),fake_dv_comb(temp_range(1):temp_range(end),2))
    end
    hold on
    plot(gauss_out{q,1});
    title(['Frame: ' num2str(q)])
    xlabel('DV position')
    ylabel('Fluorescence')
    set(gca,'FontSize',20)
    axis([-300 300 0 2500])
    saveas(gcf, ['Frame: ' num2str(q) '.tif'])
    clf
    
end