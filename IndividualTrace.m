function [Time,TotalFluo]=IndividualTrace(Transitions,Rates,Delay,MaxTime)
% [Time,TotalFluo]=IndividualTrace(Transitions,Rates,Delay,MaxTime)
%
% DESCRIPTION
% This function generates an individual trace given times for rate change
% and rates. 
%
% ARGUMENTS
% Transitions: List of time points when bursts happen. 
% Rates: List of rates during each burst (between transitions)
% Delay: Usually Delay=GeneLength/ElongationRate. Minutes for PolII to fall off after reaching
%        the first loops site.
% MaxTime: Ignore any events occuring after specified MaxTime
%
% OUTPUT
% Time: Time points for the trace
% TotalFluo: Corresponding fluorescence values
%
% Author (contact): Hernan Garcia (hggarcia@berkeley.edu)
% Created: 
% Last Updated: 


% Only need to include times where events happen
Time = [Transitions, Transitions+Delay, MaxTime];
Time = sort(Time);
% Ignore any events occuring after specified MaxTime
Time(Time>MaxTime) = [];
Fluo=zeros(length(Time),length(Transitions));

for i=1:length(Transitions)
    %Starting conditions
    for j=1:length(Time)
        %If this is not the last transition
        if i<length(Transitions)
            % Total fluorescence added (subtract later)
            if (Time(j)<=Transitions(i))
                Fluo(j,i)=0;
            elseif (Time(j)<=Transitions(i+1))
                Fluo(j,i)=Rates(i)*(Time(j)-Transitions(i));
            elseif (Time(j)<=Transitions(i+1)+Delay)
                Fluo(j,i)=Rates(i)*(Transitions(i+1)-Transitions(i));
            end
            % Subtract fluorescence
            if (Time(j)>=Transitions(i)+Delay)&&(Time(j)<=Transitions(i+1)+Delay)
                Fluo(j,i)=Fluo(j,i)-Rates(i)*(Time(j)-Transitions(i)-Delay);
            end

        %If this is the last transition - NOTE: I need to see what happens
        %if the MaxTime is within the life time of this trace
        else    
            if Time(j)<=Transitions(i)
                Fluo(j,i)=0;
            elseif (Time(j)>Transitions(i))&&(Time(j)<=Transitions(i)+Delay)
                Fluo(j,i)=Rates(i)*(Time(j)-Transitions(i));
            elseif (Time(j)>Transitions(i)+Delay)&&(Time(j)<=MaxTime)
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
