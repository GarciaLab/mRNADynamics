function checkMCPDosage(DataType)
%
%double check the offset to ensure correct mcp dosage in all data sets
close all;

if ischar(DataType)
    [allData, Prefixes, resultsFolder] = LoadMS2Sets(DataType);
else
    allData = DataType;
    DataType = inputname(1);
end

figure()
ax12 = subplot(1,3, 1)
for j = 1:length(allData)
    nc12 = allData(j).Particles.nc12;
    nc13= allData(j).Particles.nc13;
    off = allData(j).Particles.MeanOffsetVector(nc12:nc13);
    plot(1:length(off), off, 'DisplayName', ['embryo ',num2str(j+1)])
    ylabel('mean offset (au)')
    xlabel('frame')
    title('offset in nc12 for all 0DG embryos')
    hold on
    leg12 = legend(ax12);
    standardizeFigure(ax12,leg12); 
end

ax13 = subplot(1,3,2)
for j = 1:length(allData)
    nc13 = allData(j).Particles.nc13;
    nc14 = allData(j).Particles.nc14;
    off = allData(j).Particles.MeanOffsetVector(nc13:nc14);
    plot(1:length(off), off, 'DisplayName', ['embryo ',num2str(j+1)])
    ylabel('mean offset (au)')
    xlabel('frame')
    title('offset in nc13 for all 0DG embryos')
    hold on
    leg13 = legend(ax13);
    standardizeFigure(ax13,leg13); 
end

ax14 = subplot(1,3,3)
for j = 1:length(allData)
    nc14 = allData(j).Particles.nc14;
    off = allData(j).Particles.MeanOffsetVector(nc14:end);
    plot(1:length(off), off, 'DisplayName', ['embryo ',num2str(j+1)])
    ylabel('mean offset (au)')
    xlabel('frame')
    title('offset in nc14 for all 0DG embryos')
    hold on
    leg14 = legend(ax14);
    standardizeFigure(ax14,leg14); 
end

end