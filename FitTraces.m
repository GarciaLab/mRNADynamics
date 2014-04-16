
%%%%%%%%%% Fake generated data %%%%%%%

% [TimeWindow,TotalFluo]=IndividualTrace([1,2,3,20,22,25,30],[1,1,0,1,0,1,0],4,70)
% 
% x= [0:0.5:70];
% y = interp1(TimeWindow,TotalFluo,x);
% 
% yn=(0.5-rand(1,141))/5;
% 
% y=y+yn;


%%%%%%%%%%% Real Data %%%%%%%%


if exist('C:\Users\bothma\Dropbox\MS2Pausing')>0;
DropboxFolder='C:\Users\bothma\Dropbox\MS2Pausing';
else
DropboxFolder='/Users/bothma/Dropbox/MS2Pausing';
%DropboxFolder='/Users/bothma/Dropbox/MS2AnalysisJB';
end

Prefix='2014-03-01-SnaBAC_A';

Data=load([DropboxFolder,filesep,Prefix,filesep,'CompiledParticles.mat']);

ParticleNum=187;
%ParticleNum=206;
%ParticleNum=212;

y = Data.CompiledParticles(ParticleNum).Fluo/max(Data.CompiledParticles(ParticleNum).Fluo);

xx= Data.ElapsedTime(Data.CompiledParticles(ParticleNum).Frame)-...
                Data.ElapsedTime(Data.nc14);
            
            plot(xx,y)



x= [0:0.5:max(xx)];

y = 4*interp1(xx,y,x);

y(isnan(y))=0;

plot(x,y)

%     
    Transitions=[];
    Rates =[];
% 
% for i=1:length(x)
%     
%     
%     if isempty(Transitions)
%     y1=zeros(1,i);
%     else
%        [TimeWindowt,TotalFluot]=IndividualTrace(Transitions,Rates,4,70);
%        y1= interp1(TimeWindowt,TotalFluot,x);
%     end
%     
%     Lsq(i)=mean((abs(y(1:i)-y1(1:i))));
%     
%     if i>1 & (y(i)-y1(i))>0.2
%         
%         Transitions=[Transitions,x(i-1)]
%         Rates = [Rates,1]
%         
%     elseif i>1 & (y(i)-y1(i))<-0.2
%         
%         
%         Transitions=[Transitions,x(i-1)];
%         Rates = [Rates,0];
%         
%     end
%     
% end
% 
%         
% 
% plot(x,y,'r')    
% hold on, plot(x,y1,'g')    
    
Thresh=0.07;
Offset=3;

TP=1;

RR=[];

StoreTP=[];

holder=1;

while holder<=10

    
for N=TP:length(x);
    [p,S] = polyfit(TP:N, y(TP:N),1);yfit = polyval(p,TP:N);
    RR(N)=sum(abs(y(TP:N)-yfit));
end



TP=find(RR./[1:length(RR)]>=Thresh,1,'first');

if isempty(TP)
    holder=100;
else
StoreTP(holder)=TP;

holder=holder+1;
end

end


StoreTP=[1,StoreTP-Offset,length(x)];

TT=[];
m=[];
c=[];

for i=1:length(StoreTP)-1
     
    [p,S] = polyfit(StoreTP(i):StoreTP(i+1)-1,y(StoreTP(i):StoreTP(i+1)-1),1);yfit = polyval(p,StoreTP(i):StoreTP(i+1)-1);
    TT=[TT,yfit];
    m(i)=p(1);
    c(i)=p(2);
end

plot(y)
hold on
plot(TT,'r')


% 
% for N=TP:length(x);
%     [p,S] = polyfit(TP:N, y(TP:N),1);yfit = polyval(p,TP:N);
%     RR2(N)=sum(abs(y(TP:N)-yfit));
% end
% 
% TP2=find(RR2./[1:length(RR)]>=Thresh,1,'first')
% 
% for N=TP2:length(x);
%     [p,S] = polyfit(TP2:N, y(TP2:N),1);yfit = polyval(p,TP2:N);
%     RR3(N)=sum(abs(y(TP2:N)-yfit));
% end
% 
% TP3=find(RR3./[1:length(RR)]>=Thresh,1,'first')
% 
% 
% figure, plot(RR2./[1:length(RR)],'g')
% hold on
% plot(RR./[1:length(RR)],'r')
% plot(RR3./[1:length(RR3)],'m')
% 
% plot(y)
%     