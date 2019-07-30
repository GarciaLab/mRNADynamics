function dorsalSignalAnalysis(dataTypes)

d = {};
dataTypes = {'1DG', '1DG_Zeiss'};
prefixes = {};
[d{1}, prefixes{1}, resultsFolder] = LoadMS2Sets(dataTypes{1}, 'noCompiledNuclei');
[d{2}, prefixes{2}, ~] = LoadMS2Sets(dataTypes{2}, 'noCompiledNuclei');
% [d{3}, prefixes{3}, ~] = LoadMS2Sets(dataTypes{3}, 'noCompiledNuclei');

today = date;
savepath = [resultsFolder, filesep, 'figures', filesep, today];
mkdir(savepath);

ch = 1;
dynrg = cell(2, length(d));
cv = cell(2, length(d));
lens = cell(2, length(d));

for nc = 12:13
    for i = 1:length(d)
        snips = {};
        for e = 1:length(d{i})
            cp =   d{i}(e).Particles.CompiledParticles{ch};
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
                    cv{nc-11, i} = [cv{nc-11, i}, nanmedian( abs(cpp.FluoGauss3DError ./ cpp.FluoGauss3D))  ];
                    lens{nc-11, i} = [lens{nc-11, i}, length((cpp.Frame))];
                end
            end
             figure()
            montage(snips);
            saveas(gcf, [savepath, filesep, dataTypes{i}, '_nc ',num2str(nc),'_snips'])
    end
end


figure
k = 1;
axsdyn = {};
for nc = 12:13
    axsdyn{nc-11} = subplot(1, 2, nc-11);
    for i = 1:length(dynrg)
        histogram(axsdyn{nc-11}, log10(dynrg{nc-11, i}+1), 6, 'Normalization', 'count');
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
axscv = {};
k = 1;
for nc = 12:13
    axscv{nc-11} = subplot(1, 2, nc-11);
    for i = 1:length(cv)
        histogram(axscv{nc-11}, log10(cv{nc-11, i}+1), 6, 'Normalization', 'count');
        k = k + 1;
        hold on
    end
    title(num2str(nc))  
    ylabel('counts')
    xlabel('fluo / fluo error')

    legs{nc-11} = legend('1DG', '1DG_Zeiss');
    standardizeFigure(gca, []);
end


saveas(gcf, [savepath, filesep, 'cv'])



figure
axslens = {};
k = 1;
for nc = 12:13
    axslens{nc-11} = subplot(1, 2, nc-11);
    for i = 1:length(cv)
        histogram(axslens{nc-11}, lens{nc-11, i}, 'Normalization', 'count');
        k = k + 1;
        hold on
    end
    title(num2str(nc))  
    ylabel('pdf')
    xlabel('trace length')

    legs{nc-11} = legend('1DG', '1DG_Zeiss');
    standardizeFigure(gca, []);
end


saveas(gcf, [savepath, filesep, 'cv'])
   
end