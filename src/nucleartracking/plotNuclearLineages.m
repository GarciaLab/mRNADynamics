function plotNuclearLineages(schnitzcells)
%author- ar


xmax = 512;
ymax = 512;


close all;
linFig = figure();
linAx = axes(linFig);

nc = 12;
for s = 1:length(schnitzcells)
    
    schnitz = schnitzcells(s);
    
    if schnitzcells(s).cycle == nc
        frames = schnitz.frames;
        clrmp = magma(length(frames));

        for f = 1:length(frames)
            cenx = schnitz.cenx(f);
            ceny = schnitz.ceny(f);
            viscircles(linAx,[cenx,ceny], double(schnitz.len(f)),'Color', clrmp(f,:));
            %             plot(cenx, ceny, 'o', 'Color', clrmp(f,:));
            xlim([0, xmax]);
            ylim([0, ymax]);
            drawnow;
            hold on
        end
    end
    
end
    