

folder =  '/Users/bothma/Dropbox/HbMovies/2014-03-18-HbBac_A/MovieData/';

    
load([folder,'MaxNuclei.mat'],'MaxNuclei');

TotalTime=length(MaxNuclei.TimeMatrix);

FR=[];for i=1:TotalTime; FR(i) = OptimalFilterRadiusNuclearSeg(MaxNuclei.(['Time', num2str(i)]),3,8,50,0.5,20,0); i, end

FR(isnan(FR))=10;

LabelNucsCore = SegmentNucleiLiveFunction(folder,[],round(FR+2),round(8*FR));