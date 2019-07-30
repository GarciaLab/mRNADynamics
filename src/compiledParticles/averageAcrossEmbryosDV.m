function averageAcrossEmbryosDV(DataType)
close all;
if ischar(DataType)
    [allData, Prefixes, resultsFolder] = LoadMS2Sets(DataType);
else
    allData = DataType;
    DataType = inputname(1);
end
figure(4);
meanDVAll = [];

allRNAs = {[], [], []};
allDorsals = {[], [] []};

for nc = 1:3
    for e = 1:length(allData) %loop over embryos
        
        cpo = allData(e).Particles;
        cp = cpo.CompiledParticles;
        DVbinID = cpo.DVbinID;
        frames12 = cpo.cycleFrames{1};
        frames13 = cpo.cycleFrames{2};
        frames14 = cpo.cycleFrames{3};
        
        if e == 1
            frames12all = frames12;
            frames13all = frames13;
            frames14all = frames14;
        else
            if length(frames12) >= length(frames12all)
                frames12all = frames12;
            end
            if length(frames13) >= length(frames13all)
                frames13all = frames13;
            end
            if length(frames14) >= length(frames14all)
                frames14all = frames14;
            end
        end
        meanDV = cpo.meanDV{nc};
        sumDV = cpo.sumDV{nc};
        sumDV(isnan(sumDV)) = 0;
        countsDV = cpo.countsDV{nc};
        countsDV(isnan(countsDV)) = 0;
        nbins = size(meanDV, 1);
        if e == 1
            meanDVAll{nc} = meanDV;
            sumDVAll{nc} = sumDV;
            countsDVAll{nc} = countsDV;
        else
            ld = size(meanDV,2) - size(meanDVAll{nc},2);
            if ld > 0
                meanDVAll{nc} = [meanDVAll{nc}, zeros(nbins, ld)];
                sumDVAll{nc} = [sumDVAll{nc}, zeros(nbins, ld)];
                countsDVAll{nc} = [countsDVAll{nc}, zeros(nbins, ld)];
            elseif ld < 0
                meanDV = [meanDV, zeros(nbins, abs(ld))];
                sumDV = [sumDV, zeros(nbins, abs(ld))];
                countsDV = [countsDV, zeros(nbins, abs(ld))];
                
            end
            
            sumDVAll{nc} = sumDVAll{nc} + sumDV;
            countsDVAll{nc} = countsDVAll{nc} + countsDV;
            meanDVall{nc} = sumDVAll{nc}./countsDVAll{nc};
            
            
        end
        
        %      sumDV = {zeros(length(DVbinID), length(frames12)),...
        %         zeros(length(DVbinID), length(frames13)),zeros(length(DVbinID), length(frames14))};
        %     meanDV = sumDV;
        %     countsDV = {zeros(length(DVbinID), length(frames12)),...
        %         zeros(length(DVbinID), length(frames13)),zeros(length(DVbinID), length(frames14))};
        %     for nc = 12:14
        %         for p = 1:length(cp)
        %             for dv = 1:length(DVbinID)
        %                 if cp(p).cycle == nc && cp(p).dvbin == dv
        %                     sumDV{nc-11}(dv,cp(p).FramesWRTAnaphase) = sumDV{nc-11}(dv,cp(p).FramesWRTAnaphase) + cp(p).Fluo;
        %                     countsDV{nc-11}(dv,cp(p).FramesWRTAnaphase) = countsDV{nc-11}(dv,cp(p).FramesWRTAnaphase) + 1;
        %                 end
        %             end
        %         end
        %         meanDV{nc-11} = sumDV{nc-11}./countsDV{nc-11};
        %     end
        
        
        %plot by dorsal concentration
        allRNAs{nc} = [allRNAs{nc}, cpo.allRNAs{nc}];
        allDorsals{nc} = [allDorsals{nc}, cpo.allDorsals{nc}];
        
        
        
    end
    
    %%
    % dvspresent = unique(dvspresent);
    dvspresent = 1:length(DVbinID);
    cumdv = {};
    cm = magma;
    cmsize = size(cm, 1);
    cmslice = cm(1:floor(cmsize/length(dvspresent)):end, :);
    
    curAv = meanDVAll{nc};
    time = 1:size(curAv, 2);
    % figure(nc);
    for dv = 1:nbins
        cumdv{nc}(dv) = nansum(curAv(dv, :));
        subplot(3, 2, nc*2 - 1);
        if cumdv{nc}(dv)~=0
            plot(time, curAv(dv, :), 'Color', cmslice(dv,:), 'lineWidth', 3, 'DisplayName', [num2str(DVbinID(dv)*100),'% dv']);
        end
        hold on
    end
    
    title(['mean spot intensity. nc ', num2str(nc+11)]);
    xlabel('frames since anaphase');
    ylabel('intensity (au)');
    legend;
    
    subplot(3, 2, nc*2);
    cumdv{nc}(cumdv{nc}==0) = NaN;
    %     plot(1:nbins, sl, 'lineWidth', 3);
    %     cumdvnan{nc} = cumdv{nc}
    if length(dvspresent) > 2
        h = colormapline(DVbinID(dvspresent), cumdv{nc}(dvspresent),[], cmslice);
    else
        h = plot(DVbinID(dvspresent), cumdv{nc}(dvspresent), 'Color',cmslice(dvspresent, :));
    end
    %     xlim([min(dvspresent) max(dvspresent)]);
    set(h,'linewidth',3)
    title(['accumulated spot intensity. nc ', num2str(nc+11)]);
    xlabel('fraction dv');
    ylabel('intensity (au)');
    
    
    
end


figure()
for i = 1:3
    subplot(1, 3, i)
    plot(allDorsals{i}, allRNAs{i}, 'o');
    set(gca, 'YScale', 'log');
    xlabel('dorsal concentration (au)');
    ylabel('accumulated RNA per particle');
    title(['nc',num2str(i+11)]);
    standardizeFigure(gca, []);
end


end

