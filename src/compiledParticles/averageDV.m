function averageDV(resultsFolder)

load([resultsFolder, filesep, 'CompiledParticles.mat']);
load([resultsFolder, filesep, 'APDivision.mat']);

% meanDV12 = zeros(length(DVbinID));

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
    frames14 = 1:(length(ElapsedTime)-min14);
    
    sumDV = {zeros(length(DVbinID), length(frames11)), zeros(length(DVbinID), length(frames12)),...
        zeros(length(DVbinID), length(frames13)),zeros(length(DVbinID), length(frames14))};
    meanDV = sumDV;
    countsDV = {ones(length(DVbinID), length(frames11)), ones(length(DVbinID), length(frames12)),...
        ones(length(DVbinID), length(frames13)),ones(length(DVbinID), length(frames14))};
    for nc = 12:14
        for p = 1:length(cp)
            for dv = 1:length(DVbinID)
                if cp(p).cycle == nc && cp(p).dvbin == dv
                    sumDV{nc-10}(dv,cp(p).FramesWRTAnaphase) = sumDV{nc-10}(dv,cp(p).FramesWRTAnaphase) + cp(p).Fluo;
                    countsDV{nc-10}(dv,cp(p).FramesWRTAnaphase) = countsDV{nc-10}(dv,cp(p).FramesWRTAnaphase) + 1;
                end
            end
        end
        meanDV{nc-10} = sumDV{nc-10}./countsDV{nc-10};
    end
end

save([resultsFolder, filesep, 'CompiledParticles.mat'],'meanDV', 'countsDV', 'sumDV','-append');


end