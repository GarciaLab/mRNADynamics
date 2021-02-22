

folder =  '/Users/bothma/Dropbox/HbMovies/2014-03-18-HbBac_A/MovieData/';

    
load([folder,'MaxNuclei.mat'],'MaxNuclei');

load([folder,'LabelNucsCore.mat'],'LabelNucsCore');


TotalTime = length(fieldnames(MaxNuclei))-1;
 
%%%%%%%%%%% NucDilate %%%%%%%%%%%%%%

    str='';

if exist([folder,'LabNucDilate.mat'])>0
    
    str = input('LabNucDilate already exists, overwrite? (y/n):','s');
    
end
    
  if strcmp('y',str)|~(exist([folder,'LabNucDilate.mat']))
    
    for i=1:TotalTime

        [lengthh,widthh]=size(MaxNuclei.Time1);
        
      LabNucDilate.(['Time', num2str(i)]).Image  = segmentnucleiExpandLive(MaxNuclei.(['Time', num2str(i)]),LabelNucsCore.(['Time', num2str(i)]).ImageZ,50,20,ones(lengthh,widthh),50,20,10);

    end
    
          save([folder,'LabNucDilate.mat'],'LabNucDilate');
        
    else
        disp('LabNucDilate already exists, loading file')
        
        load([folder,'LabNucDilate.mat'],'LabNucDilate');
    end
    
    
    NumNucs=[];

for i=1:TotalTime
    
    NumNucs=[unique(LabelNucsCore.(['Time', num2str(i)]).ImageZ(:)); NumNucs];
end


NumNucs=max(NumNucs);

ColorMa=colormap(jet(50+NumNucs));

ColorMaRand=ColorMa(randperm(50+NumNucs),:);

for i=1:TotalTime
    
    CM=label2rgb(LabNucDilate.(['Time', num2str(i)]).Image,ColorMaRand,[1,1,1]);
    
DI = cast(bsxfun(@times,double(CM)/(255),double(imadjust(MaxNuclei.(['Time', num2str(i)])))),class(MaxNuclei.(['Time', num2str(i)])));
    
imshow(DI);
 pause(0.1)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






% TotalTime=length(MaxNuclei.TimeMatrix);
%
% FR=[];
% for i=1:TotalTime;
%     FR(i) = OptimalFilterRadiusNuclearSeg(MaxNuclei.(['Time', num2str(i)]),3,15,100,0.5,20,0);
%     i,
% end
% 
% FR(isnan(FR))=10;
%
% LabelNucsCore = SegmentNucleiLiveFunction(folder,[],round(FR+2),round(8*FR));