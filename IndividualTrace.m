function [TimeWindow,TotalFluo]=IndividualTrace(Transitions,Rates,Delay,MaxTime)

%This function generates an individual trace given times for rate change
%and rates

TimeWindow=linspace(0,MaxTime);
Fluo=zeros(length(TimeWindow),length(Transitions));

for i=1:length(Transitions)
       
    %Starting conditions
    for j=1:length(TimeWindow)
        %If this is not the last transition
        if i<length(Transitions)

            if TimeWindow(j)<=Transitions(i)
                Fluo(j,i)=0;
            elseif (TimeWindow(j)>Transitions(i))&(TimeWindow(j)<=Transitions(i)+Delay)
                Fluo(j,i)=Rates(i)*(TimeWindow(j)-Transitions(i));

            elseif (TimeWindow(j)>Transitions(i)+Delay)&(TimeWindow(j)<=Transitions(i+1))
                Fluo(j,i)=Rates(i)*Delay;

            elseif (TimeWindow(j)>Transitions(i+1))&(TimeWindow(j)<=Transitions(i+1)+Delay)
                Fluo(j,i)=Rates(i)*(-TimeWindow(j)+Delay+Transitions(i+1));

            end
        %If this is the last transition - NOTE: I need to see what happens
        %if the MaxTime is within the life time of this trace
        else    
            if TimeWindow(j)<=Transitions(i)
                Fluo(j,i)=0;
            elseif (TimeWindow(j)>Transitions(i))&(TimeWindow(j)<=Transitions(i)+Delay)
                Fluo(j,i)=Rates(i)*(TimeWindow(j)-Transitions(i));

            elseif (TimeWindow(j)>Transitions(i)+Delay)&(TimeWindow(j)<=MaxTime)
                Fluo(j,i)=Rates(i)*Delay;

%             elseif (TimeWindow(j)>Transitions(i+1))&(TimeWindow(j)<=Transitions(i+1)+Delay)
%                 Fluo(i,j)=Rates(i)*(-TimeWindow(j)+Delay+Transitions(i+1));

            end
        end
        
      
    end
end
    
%Finally, truncate everything that is below zero
Fluo(Fluo<0)=0;
        


TotalFluo=sum(Fluo,2);

% figure(2)
% hold all
% for i=1:size(Fluo,2)
%     plot(TimeWindow,Fluo(:,i),'-')
% end
% hold off
