function chi2=lsqnonlinFitIndividualFluorescenceCurve(TimeData,FluoData,Delay,nSteps,x0)

%Gives the chi square of the data to the fit form in order to do a fit with
%lsqnonlin

%Starting conditions
%nSteps=x0(1);
Transitions0=x0(1:nSteps);
Rates0=x0(1+nSteps:end);

%Get the predicted shape
[TimePrediction,FluoPrediction]=IndividualTrace(Transitions0,Rates0,Delay,max(TimeData));
FluoPrediction=interp1(TimePrediction,FluoPrediction,TimeData);

% Fluo=zeros(length(TimeData),length(Transitions0));
% 
% for i=1:length(Transitions0)
%        
%     %Starting conditions
%     for j=1:length(TimeData)
%         %If this is not the last transition
%         if i<length(Transitions0)
% 
%             if TimeData(j)<=Transitions0(i)
%                 Fluo(i,j)=0;
%             elseif (TimeData(j)>Transitions0(i))&(TimeData(j)<=Transitions0(i)+Delay)
%                 Fluo(i,j)=Rates0(i)*(TimeData(j)-Transitions0(i));
% 
%             elseif (TimeData(j)>Transitions0(i)+Delay)&(TimeData(j)<=Transitions0(i+1))
%                 Fluo(i,j)=Rates0(i)*Delay;
% 
%             elseif (TimeData(j)>Transitions0(i+1))&(TimeData(j)<=Transitions0(i+1)+Delay)
%                 Fluo(i,j)=Rates0(i)*(-TimeData(j)+Delay+Transitions0(i+1));
% 
%             end
%         %If this is the last transition - NOTE: I need to see what happens
%         %if the MaxTime is within the life time of this trace
%         else    
%             if TimeData(j)<=Transitions0(i)
%                 Fluo(i,j)=0;
%             elseif (TimeData(j)>Transitions0(i))&(TimeData(j)<=Transitions0(i)+Delay)
%                 Fluo(i,j)=Rates0(i)*(TimeData(j)-Transitions0(i));
% 
%             elseif (TimeData(j)>Transitions0(i)+Delay)&(TimeData(j)<=max(TimeData))
%                 Fluo(i,j)=Rates0(i)*Delay;
% 
%             end
%         end
%         
%        
%       
%     end
% end
% 
% %Truncate anything smaller than zero
% Fluo(Fluo<0)=0;
%             
% 
% FluoPrediction=sum(Fluo);





%Calculate chi^2
try
    chi2=(FluoData-FluoPrediction).^2;
catch
    1+1;
end

