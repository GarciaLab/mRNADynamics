function dorsalSignalAnalysis(dataTypes)

d = {};
dataTypes = {'0DG', '1DG', '1DG_Zeiss'};
prefixes = {};
[d{1}, prefixes{1}, resultsFolder] = LoadMS2Sets(dataTypes{1});
[d{2}, prefixes{2}, ~] = LoadMS2Sets(dataTypes{2});
[d{3}, prefixes{3}, ~] = LoadMS2Sets(dataTypes{3});

today = date;
savepath = [resultsFolder, filesep, 'figures', filesep, today];
mkdir(savepath);

ch = 1;
dynrg = cell(2, length(d));
cv = cell(2, length(d));


for nc = 12:13
    for i = 1:length(d)
        snips = {};
        for e = 1:length(d{i})
            cp =   d{i}(e).Particles.CompiledParticles{ch}
            sc = d{i}(e).Particles.schnitzcells;
                snipdirpath = [resultsFolder, filesep,prefixes{i}{e},filesep, 'ParticleTraces'];
                snipDir = dir([snipdirpath, filesep, '*nc',num2str(nc),'*.tif']);
                for j = 1:length(snipDir)
                    snips = [snips, imread([snipdirpath, filesep, snipDir(j).name])];
                end
            end
            
            for p = 1:length(cp)
                cpp = cp(p);
                if cpp.cycle == nc & cpp.Approved & sc(cpp.schnitz).Approved
                    dynrg{nc-11, i} = [dynrg{nc-11, i}, max(cpp.FluoGauss3D) / abs(min(cpp.FluoGauss3D))];
                    cv{nc-11, i} = [cv{nc-11, i}, nanmean( abs(cpp.FluoGauss3DError ./ cpp.FluoGauss3D))  ];
                end
            end
             figure()
            montage(snips);
            saveas(gcf, [savepath, filesep, dataTypes{i}, '_nc ',num2str(nc),'_snips'])
    end
end


figure
% axs = {};
k = 1;
for nc = 12:13
    axs{nc-11} = subplot(1, 2, nc-11);
    for i = 2:length(dynrg)
        histogram(axs{nc-11}, log10(dynrg{nc-11, i}+1), 6, 'Normalization', 'pdf');
        k = k + 1;
        hold on
    end
    title(num2str(nc))  
    ylabel('counts')
    xlabel('log dynamic range')

    legs{nc-11} = legend('1DG', '1DG_Zeiss');
    standardizeFigure(gca, []);
end


saveas(gcf, [savepath, filesep, 'dynamicRange'])



figure
% axs = {};
k = 1;
for nc = 12:13
    axs{nc-11} = subplot(1, 2, nc-11);
    for i = 2:length(cv)
%             subplot(2, 3, k)

        histogram(axs{nc-11}, log10(cv{nc-11, i}+1), 6, 'Normalization', 'pdf');
        k = k + 1;
        hold on
%     title([dataTypes{i}, ' ',num2str(nc)])
%     xlim([.1, 1.5])
%     ylim([.1, 5])
%     legs{nc-11} = legend('1DG', '1DG_Zeiss');
%     standardizeFigure(gca, []);
    end
    title(num2str(nc))  
    ylabel('counts')
    xlabel('fluo / fluo error')

    legs{nc-11} = legend('1DG', '1DG_Zeiss');
    standardizeFigure(gca, []);
end


saveas(gcf, [savepath, filesep, 'cv'])
   
end