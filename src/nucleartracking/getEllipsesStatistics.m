function out = getEllipsesStatistics(Ellipses)

cenxList = [];
cenyList = [];
semiMajList = [];
semiMinList = [];
orientationAngleList = [];


for f = 1:length(Ellipses)
    
    
    cenxList = [cenxList; Ellipses{f}(:, 1)];
    cenyList = [cenyList; Ellipses{f}(:, 2)];
    semiMajList = [semiMajList; Ellipses{f}(:, 3)];
    semiMinList = [semiMinList; Ellipses{f}(:, 4)];
    orientationAngleList = [orientationAngleList; Ellipses{f}(:, 5)];
    
    
end


out = table(cenxList, cenyList, semiMajList, semiMinList, orientationAngleList); 