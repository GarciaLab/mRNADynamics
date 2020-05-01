function [schnitzcells, Ellipses, CompiledParticles, Particles] =...
    ...
    breakUpSchnitzesAtMitoses(schnitzcells, Ellipses, ncs, nFrames, varargin)

p = false;
cp = false;
Particles = [];
CompiledParticles = [];

if ~isempty(varargin)
    if length(varargin) == 1
        Particles = varargin{1};
        p = true;
    elseif length(varargin) == 2
        Particles = varargin{1};
        CompiledParticles = varargin{2};
        p = true;
        cp = true;
    end
end

cycleFrames = nan(1,nFrames);

for i = 1:length(ncs)
    if i==14
        cycleFrames(ncs(i):end) = 14;
    elseif ncs(i) == 0
        cycleFrames(1:ncs(i+1)) = i;
    else
        cycleFrames(ncs(i):ncs(i+1)) = i;
    end
end


tempSchnitzcells = schnitzcells;
nNuclei = length(schnitzcells);

j = 1;
for s = 1:nNuclei
    
    schnitzcells(s).deleteMe = false;
    sc  = schnitzcells(s);
    tempSchnitzcells(s).deleteMe = false;
    schnitzcells(s).deleteMe = false;
    
    if p
        pInd = find(Particles.Nucleus == s);
    end
    if cp
        cpInd = find(CompiledParticles.Nucleus == s);
    end
    
    hasncs = unique(cycleFrames(sc.frames));
    
    if length(hasncs) > 1
        
        %         keyboard
        for i = 1:length(hasncs)
            newInd = nNuclei + j;
            tempSchnitzcells(newInd) = sc;
            newFrames = cycleFrames(sc.frames) == hasncs(i);
            
            tempSchnitzcells(newInd).frames = sc.frames(newFrames);
            tempSchnitzcells(newInd).cenx = sc.cenx(newFrames);
            tempSchnitzcells(newInd).ceny = sc.ceny(newFrames);
            tempSchnitzcells(newInd).cellno = sc.cellno(newFrames);
            tempSchnitzcells(newInd).deleteMe = false;
            
            if isfield(tempSchnitzcells, 'len')
                tempSchnitzcells(newInd).len = sc.len(newFrames);
            end
            
            if isfield(tempSchnitzcells, 'Fluo')
                tempSchnitzcells(newInd).Fluo = sc.Fluo(newFrames, :);
            end
            tempSchnitzcells(newInd).cellno = sc.cellno(newFrames);
            if isfield(tempSchnitzcells, 'APpos')
                tempSchnitzcells(newInd).APpos = sc.APpos(newFrames);
                
            end
            
            if isfield(tempSchnitzcells, 'DVpos')
                tempSchnitzcells(newInd).DVpos = sc.DVpos(newFrames);
            end
            if isfield(tempSchnitzcells, 'FrameApproved')
                tempSchnitzcells(newInd).FrameApproved = sc.FrameApproved(newFrames);
            end
            
            if isfield(tempSchnitzcells, 'FluoTimeTrace')
                
                tempSchnitzcells(newInd).FluoTimeTrace = sc.FluoTimeTrace(newFrames);
            end
            
            if cp
                CompiledParticles(cpInd).Nucleus = newInd;
                CompiledParticles(cpInd).schnitz = newInd;
            end
            if p
                Particles(pInd).Nucleus= newInd;
            end
            
            j = j+1;
        end
        tempSchnitzcells(s).deleteMe = true;
        schnitzcells(s).deleteMe = true;
        
    end
    
end

tempSchnitzcells([tempSchnitzcells.deleteMe]) = [];
if isfield(tempSchnitzcells, 'Valid')
    tempSchnitzcells = rmfield(tempSchnitzcells, {'deleteMe', 'StitchedTo', 'StitchedFrom', 'Valid'});
else
    tempSchnitzcells = rmfield(tempSchnitzcells, {'deleteMe'});
end

schnitzcells = tempSchnitzcells;

Ellipses = addSchnitzIndexToEllipses(Ellipses, schnitzcells);

