function binDorsal(DataType, varargin)

cpFlag = false;
if ~isempty(varargin)
    cpFlag = true;
end


if ~cpFlag
    [allData, Prefixes, resultsFolder] = LoadMS2Sets(DataType, 'noCompiledNuclei');
else
    [~, resultsFolder, Prefixes] = getDorsalPrefixes(DataType);
end

mind = 1;
maxd = 1500;
nBins = 20;

% dlfluobins = mind:round(maxd/mind):maxd; %this is appropriate for taking instantaneous dorsal at ~50% through nc12 on the sp8
dlfluobins = 0:250:4500;
% dlfluobins = [0:20:500, 501:250:3000];
% dlfluobins = [0:10:249,250:50:499,500:250:4000];
% dlfluobins = logspace(0, 3.6, 20);
%% the fancy method to create pseudo dv spatial bins

fspace = @(f, f1, mi, ma, n)  f1( linspace( f(mi), f(ma), n ) ); %generalization of the logspace and linspace functions
sigma = 1/(maxd*sqrt(2*pi));
c = 1/(2*sigma^2);
p = 2;
% dorsalGradient = @(x) maxd*exp(-c.*x.^p );
% inverseGradient = @(x) ((-1/c).* log(x./maxd)).^(1/p);

calib = 1500/1080;
B = 520*calib; %constant offset
M = -89*calib; %linear offset
B = 0;
M = 0;
c = 5; %related to sigma
dorsalGradient = @(x) (maxd*exp(-c*x.^p ) + B + M.*abs(x));
inverseGradient = @(x) (-(1/c).* log( (x-B-M.*abs(x)) /maxd)).^(1/p);


uDorsalGradient = @(x) -(maxd*exp(-c*x.^p ) + B + M.*abs(x)) + maxd;
uInverseGradient = @(x) (-(1/c).* log( -(x - maxd - B - M.*abs(x)) /maxd)).^(1/p);


dv = linspace(0,1, 20);
% figure
% plot(dv, dorsalGradient(dv))
maxd2 = max(uDorsalGradient(dv));
% dlfluobins = fspace(uInverseGradient,uDorsalGradient, mind, maxd2, nBins);
% dlfluobins = [linspace(71,470, round(nBins/2)), linspace(471, 3027, round(nBins/2))];

% dlfluobins = [71:10:470:, linspace(471, 3027, round(nBins/2))];

% figure
% bar(dlfluobins)
% figure()
save([resultsFolder,filesep,DataType,filesep,'dlfluobins.mat'], 'dlfluobins');

allDorsal = [];
dlfluobincounts = zeros(1, length(dlfluobins));

for e = 1:length(Prefixes)
    if ~cpFlag
        schnitzcells = allData(e).Particles.schnitzcells;
        CompiledParticles = allData(e).Particles.CompiledParticles;
        approvedIndices = [schnitzcells.Approved];
        
        for s = find(approvedIndices)
            dif = schnitzcells(s).FluoFeature - dlfluobins;
            [~,dlfluobin] = min(dif(dif>0));
            if ~isempty(dlfluobin)
                schnitzcells(s).dorsalFluoBin = single(dlfluobin);
                dlfluobincounts(dlfluobin) = dlfluobincounts(dlfluobin) + 1;
            else
                schnitzcells(s).dorsalFluoBin = NaN;
            end
        end
        
        
        save([resultsFolder,filesep,Prefixes{e},filesep,Prefixes{e},'_lin.mat'], 'schnitzcells')
        
        ch=1;
        for p = 1:length(CompiledParticles{ch})
            schnitzInd = CompiledParticles{ch}(p).schnitz;
            if schnitzcells(schnitzInd).Approved
                CompiledParticles{ch}(p).dorsalFluoBin = single(schnitzcells(schnitzInd).dorsalFluoBin);
            end
        end
        
        save([resultsFolder,filesep,Prefixes{e},filesep,'CompiledParticles.mat'], 'CompiledParticles', '-append');
        
    else
        load([resultsFolder,filesep,Prefixes{e},filesep,'compiledProject.mat'], 'compiledProject');
        for s = 1:length(compiledProject)
            dlfluo = compiledProject(s).dorsalFluoFeature;
            if ~isnan(dlfluo)
                allDorsal = [allDorsal, dlfluo];
            end
            if ~strcmpi(DataType, '1Dg')
                dif = dlfluo - dlfluobins;
            else
                dif = (2*dlfluo) - dlfluobins;
            end
            [~,dlfluobin] = min(dif(dif>0));
            if ~isempty(dlfluobin)
                compiledProject(s).dorsalFluoBin = single(dlfluobin);
                dlfluobincounts(dlfluobin) = dlfluobincounts(dlfluobin) + 1;
            else
                compiledProject(s).dorsalFluoBin = NaN;
            end
        end
        
        save([resultsFolder,filesep,Prefixes{e},filesep,'compiledProject.mat'], 'compiledProject');
        
    end
end
%
% sortedDorsal = sort(allDorsal);
% figure
% bar(sortedDorsal)
% xlabel('nucleus');
% ylabel('dl intensity');
% cumDorsal = cumtrapz(sortedDorsal);
% figure
% bar(cumDorsal);
% figure
% diffDorsal = diff(sortedDorsal);
% diff2Dorsal = diff(diffDorsal);
% bar(diffDorsal);
% figure
% bar(diff2Dorsal);
% corner = 115;
% firstpiece = linspace(sortedDorsal(1),sortedDorsal(corner),  round(nBins/2));
% secondpiece = linspace(sortedDorsal(corner)+1,sortedDorsal(end),  round(nBins/2));
% dlfluobins = [firstpiece, secondpiece]
% dlfluobins = [linspace(71,160, round(nBins/2)), linspace(161, 3027, round(nBins/2)]


save([resultsFolder,filesep,DataType,filesep,'dlfluobincounts.mat'], 'dlfluobincounts');

end