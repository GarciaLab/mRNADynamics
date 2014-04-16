
%%%%%%%%%% Fake generated data %%%%%%%

[TimeWindow,TotalFluo]=IndividualTrace([1,2,3,20,22,25,30,35],[1,0,1,0,1,0,1,0],4,70)

x= [0:0.5:70];
y = interp1(TimeWindow,TotalFluo,x);

yn=(0.5-rand(1,141))/10;

y=y+yn;


%%%%%%%%%%% Real Data %%%%%%%%

% 
% if exist('C:\Users\bothma\Dropbox\MS2Pausing')>0;
% DropboxFolder='C:\Users\bothma\Dropbox\MS2Pausing';
% else
% DropboxFolder='/Users/bothma/Dropbox/MS2Pausing';
% %DropboxFolder='/Users/bothma/Dropbox/MS2AnalysisJB';
% end
% 
% Prefix='2014-03-01-SnaBAC_A';
% 
% Data=load([DropboxFolder,filesep,Prefix,filesep,'CompiledParticles.mat']);
% 
% y = Data.CompiledParticles(187).Fluo/max(Data.CompiledParticles(187).Fluo);
% 
% xx= Data.ElapsedTime(Data.CompiledParticles(187).Frame)-...
%                 Data.ElapsedTime(Data.nc14);
%             
%             plot(xx,y)
% 
% 
% 
% x= [0:0.5:max(xx)];
% 
% y = interp1(xx,y,x);
% 
% 
% 
% plot(x,y)
% 
%     
%     Transitions=[];
%     Rates =[];

for i=1:141
    
    
    if isempty(Transitions)
    y1=zeros(1,i);
    else
       [TimeWindowt,TotalFluot]=IndividualTrace(Transitions,Rates,4,70);
       y1= interp1(TimeWindowt,TotalFluot,x);
    end
    
    Lsq(i)=mean((abs(y(1:i)-y1(1:i))));
    
    if i>1 & (y(i)-y1(i))>0.2
        
        Transitions=[Transitions,x(i-1)]
        Rates = [Rates,1]
        
    elseif i>1 & (y(i)-y1(i))<-0.2
        
        
        Transitions=[Transitions,x(i-1)];
        Rates = [Rates,0];
        
    end
    
end

        

plot(x,y,'r')    
hold on, plot(x,y1,'g')    
    
    
    
    