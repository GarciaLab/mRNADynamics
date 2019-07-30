function averageDV(Prefix)


[~,~,DropboxFolder,~, PreProcPath,...
    ~, ~, ~, ~, ~,~] = readMovieDatabase(Prefix);
resultsFolder = [DropboxFolder, filesep, Prefix];

load([resultsFolder, filesep, 'CompiledParticles.mat']);
load([resultsFolder, filesep, 'APDivision.mat']);

dvspresent = [];
for ch = 1:length(CompiledParticles)
    cp = CompiledParticles{ch};
    %make common cycle
    ap11 = APDivision(11,:);
    ap12 = APDivision(12,:);
    ap13 = APDivision(13,:);
    ap14 = APDivision(14,:);
    min11 = min(ap11(ap11~=0));
    max11 = max(ap11(ap11~=0));
    min12 = min(ap12(ap12~=0));
    max12 = max(ap12(ap12~=0));
    min13 = min(ap13(ap13~=0));
    max13 = max(ap13(ap13~=0));
    min14 = min(ap14(ap14~=0));
    frames11 = 1:min12-min11;
    frames12 = 1:min13-min12;
    frames13 = 1:min14-min13;
    frames14 = 1:(length(ElapsedTime)-min14) ;
    cycleFrames = {frames12, frames13, frames14};
    
    sumDV = {zeros(length(DVbinID), length(frames12)),...
        zeros(length(DVbinID), length(frames13)),zeros(length(DVbinID), length(frames14))};
    meanDV = sumDV;
    countsDV = {zeros(length(DVbinID), length(frames12)),...
        zeros(length(DVbinID), length(frames13)),zeros(length(DVbinID), length(frames14))};
    
    
    for nc = 12:14
        for p = 1:length(cp)
            for dv = 1:length(DVbinID)
                if cp(p).cycle == nc && cp(p).dvbin == dv
                    dvspresent = [dvspresent, dv];
                    try
                    sumDV{nc-11}(dv,cp(p).FramesWRTAnaphase) = sumDV{nc-11}(dv,cp(p).FramesWRTAnaphase) + cp(p).Fluo;
                    countsDV{nc-11}(dv,cp(p).FramesWRTAnaphase) = countsDV{nc-11}(dv,cp(p).FramesWRTAnaphase) + 1;
                    catch
                        %if a particle exists in two cycles, just skip
                        %it for simplicity. 
                    end
                end
            end
        end
        meanDV{nc-11} = sumDV{nc-11}./countsDV{nc-11};
    end
end

dvspresent = unique(dvspresent);

cumdv = {};
close all;
cm = magma(length(dvspresent));
for c = 1:3
    curAv = meanDV{c};
    nbins = size(curAv, 1);
    time = 1:size(curAv, 2);
    figure(c);
    for dv = 1:nbins
        cumdv{c}(dv) = nansum(curAv(dv, :));
        subplot(1, 2, 1);
        if cumdv{c}(dv)~=0
            plot(time, curAv(dv, :), 'Color', cm(dv,:), 'lineWidth', 3, 'DisplayName', [num2str(DVbinID(dv)*100),'% dv']);
        end
        hold on
    end
    
    title(['mean spot intensity. nc ', num2str(c+11)]);
    xlabel('frames since anaphase');
    ylabel('intensity (au)');
    legend;
    
    subplot(1, 2, 2);
%     plot(1:nbins, sl, 'lineWidth', 3);
%     cumdvnan{c} = cumdv{c}
    if length(dvspresent) > 2
        h = colormapline(DVbinID(dvspresent), cumdv{c}(dvspresent),[], cm);
    else
        h = plot(DVbinID(dvspresent), cumdv{c}(dvspresent), 'Color',cm(1, :));   
    end
%     xlim([min(dvspresent) max(dvspresent)]);
    set(h,'linewidth',3) 
    title(['accumulated spot intensity. nc ', num2str(c+11)]);
    xlabel('fraction dv');
    ylabel('intensity (au)');    
    
end



save([resultsFolder, filesep, 'CompiledParticles.mat'],'meanDV', 'countsDV', 'sumDV','cumdv','cycleFrames','-append');


end