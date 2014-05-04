function [XOut,FOut,Chi2] = lsqlinAutoFitFluorescenceCurve(F,T,XPrevious,Delay,TimeRes,FErr)

% F : Fluorescence values to fit
% T : Allowed indices to vary rate
% Memory : Number of indices
% XPrevious : Rates at previous times which will contribute to fluorescence
% FErr : Fluorescence error

%%% DEFAULT INPUTS
if nargin < 3
    XPrevious=0;
end
if nargin < 4
    Delay=4;
end
if nargin < 5
    TimeRes=0.5;
end
if nargin < 6
    FErr=1;
end
Memory=floor(Delay/TimeRes);
ExcessF=(Delay-(Memory*TimeRes))/TimeRes;
% ExcessF=0;

%%% CONDITION INPUTS
% Check that no rates are to be set where we don't have data
if ~isempty(T) && T(end)>length(F)
    disp('ERROR: Asked to set a rate at a time without fluorescence data')
end
% If XPrevious empty, assume 0
if isempty(XPrevious)
    XPrevious=0;
end
% Fix dimensions
if size(T,2)>1
    T=T';
end
if size(F,2)>1
    F=F';
end
if size(XPrevious,2)>1
    XPrevious=XPrevious';
end
% XPrevious
% If first T is greater than 1, assume the rates before are equal to the
% last previous rate, and don't fit fluos there
% NumToAdd=T(1)-1;
if isempty(T)
    NumToAdd=length(F);
else
    NumToAdd=T(1)-1;
end

XPreviousAug=XPrevious;
XPreviousAug(end+1:end+NumToAdd,1)=XPrevious(end);

%%% FIT CALCULATIONS
NF=length(F);
NT=length(T);
NPre=length(XPrevious);
NPreAug=NPre+NumToAdd;

%%% PREVIOUS CONTRIBUTION
APre=ElongationMatrix(NPreAug+Memory,NPreAug,Memory,ExcessF);
% Generate previous fluorescence from augmented XPrevious
FPre=APre*XPreviousAug;
% Cut off values of FPre that aren't in window
FPre=FPre(NPre+1:end);
% This amount of fluorescence already accouted for, so subtract
FFit=F;
MaxIdx=min(length(FPre),length(FFit));
FFit(1:MaxIdx)=FFit(1:MaxIdx)-FPre(1:MaxIdx);

%%% FITTING
% Create elongation matrix
A=ElongationMatrix(NF,NF,Memory,ExcessF);
% For fitting, cut off fluorescence before first free rate...
FFit=FFit(1+NumToAdd:end,1);
% AFitStep1=A(1+NumToAdd:end,1+NumToAdd:end);
AFitStep1=A(1+NumToAdd:end,:);
% ...and combine columns for constrained rates
AFit=[];
for i=1:NT
    MinCol=T(i);
    if i==NT
        MaxCol=NF;
    else
        MaxCol=T(i+1)-1;
    end
    AFit(:,i)=sum(AFitStep1(:,MinCol:MaxCol),2);
end

% Do the fit
[XFitted,resnorm,residual]=lsqnonneg(AFit,FFit);
% size(XFitted)
% if isempty(AFit)
%     XFitted=[];
% else
%     options=optimoptions(@lsqlin,'Algorithm','active-set');
% %     size(FFit)
% %     size(AFit)
%     n=size(AFit,2);
%     lb=zeros(n,1);
%     [XFitted,resnorm,residual]=lsqlin(AFit,FFit,[],[],[],[],lb,[],[]);
% %     [XFitted,resnorm,residual]=lsqlin(AFit,FFit,-eye(n,n),zeros(n,1),[],[],[],[],[],options);
% end
FFitted=AFit*XFitted;

%%% OUTPUT FORMAT
% Pad rates with initial rate at start
XOut=zeros(NF,1);
if NumToAdd>0
    XOut(1:NumToAdd)=XPrevious(end);
end
% Fill in other values of XOut
for i=1:NT
    MinIdx=T(i);
    if i==NT
        MaxIdx=NF;
    else
        MaxIdx=T(i+1)-1;
    end
    XOut(MinIdx:MaxIdx)=XFitted(i);
end

% Pad fluo with zeros
FOut=zeros(NF,1);
FOut(1+NumToAdd:end)=FFitted;
% Add in previous contribution
MaxIdx=min(length(FPre),length(FOut));
FOut(1:MaxIdx)=FOut(1:MaxIdx)+FPre(1:MaxIdx);

% Goodness of fit
Chi2=sum((FOut-F).^2)./FErr.^2;

end