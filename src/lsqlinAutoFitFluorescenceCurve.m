function [XOut,FOut,Chi2] = lsqlinAutoFitFluorescenceCurve(F,T,XPrevious,Delay,TimeRes,FErr)
% lsqlinAutoFitFluorescenceCurve: Helper function for
% AutoFitFluorescenceCurveLinear
% Finds best fit using linear least squares, given fluorescence data and
% a set of free rates (rates to vary)
%
% F : Fluorescence values to fit
% T : Indices of free rates (referenced to start of window)
% XPrevious : Rates at times before window but within memory
% Delay : Elongation time (in minutes, not indices)
% TimeRes: Time between points
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
% Memory : elongation time (in indices, rounded down)
Memory=floor(Delay/TimeRes);
% ExcessF : amount of fluorescence leaking into next time point after
% Memory, as a fraction of full amount
ExcessF=(Delay-(Memory*TimeRes))/TimeRes;

%%% CONDITION INPUTS
% Check that we have data at all times to fit
if ~isempty(T) && T(end)>length(F)
    disp('ERROR: Can only fit rates at times with fluo data')
end
% If XPrevious empty, assume 0
if isempty(XPrevious)
    XPrevious=0;
end
% Force column vectors
if size(T,2)>1
    T=T';
end
if size(F,2)>1
    F=F';
end
if size(XPrevious,2)>1
    XPrevious=XPrevious';
end

%%% PREVIOUS RATES
% Idea: if T(1) > 1, set rates for times [1,T(1)] to XPrevious(end)
% If T empty, set all rates to XPrevious(end)
if isempty(T)
    NumToAdd=length(F);
else
    NumToAdd=T(1)-1;
end
% XPreviousAug : XPrevious, augmented with rates for times [1,T(1)]
XPreviousAug=XPrevious;
XPreviousAug(end+1:end+NumToAdd,1)=XPrevious(end);

%%% FIT CALCULATIONS
% NF : number of fluorescence values given
NF=length(F);
% NT : number of transitions
NT=length(T);
% NPre : number of rates given before window
NPre=length(XPrevious);
% NPreAug : number of rates before T(1)
NPreAug=NPre+NumToAdd;

%%% PREVIOUS CONTRIBUTION
% Generate previous fluorescence from XPreviousAug
APre=ElongationMatrix(NPreAug+Memory,NPreAug,Memory,ExcessF);
FPre=APre*XPreviousAug;
% Cut off values of FPre that aren't in window
FPre=FPre(NPre+1:end);
% Previous fluorescence already accouted for, so subtract from F before
% fitting
FFit=F;
MaxIdx=min(length(FPre),length(FFit));
FFit(1:MaxIdx)=FFit(1:MaxIdx)-FPre(1:MaxIdx);

%%% FITTING
% Create elongation matrix
A=ElongationMatrix(NF,NF,Memory,ExcessF);
% For fitting, cut off fluorescence before T(1)
FFit=FFit(1+NumToAdd:end,1);
% Step 1 because we will process more
AFitStep1=A(1+NumToAdd:end,:);

% Combine columns for constrained rates
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
FFitted=AFit*XFitted;

% TEST CODE: trying to do constrained lsqlin fit (upper bound on rate)
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
% FFitted=AFit*XFitted;

%%% PUT INTO OUTPUT FORMAT
% Rates
% Previous Rates
XOut=zeros(NF,1);
if NumToAdd>0
    XOut(1:NumToAdd)=XPrevious(end);
end
% Fitted Rates
for i=1:NT
    MinIdx=T(i);
    if i==NT
        MaxIdx=NF;
    else
        MaxIdx=T(i+1)-1;
    end
    XOut(MinIdx:MaxIdx)=XFitted(i);
end

% Fluo
% Fitted Fluo
FOut=zeros(NF,1);
FOut(1+NumToAdd:end)=FFitted;
% Previous Fluo
MaxIdx=min(length(FPre),length(FOut));
FOut(1:MaxIdx)=FOut(1:MaxIdx)+FPre(1:MaxIdx);

% Goodness of fit
Chi2=sum((FOut-F).^2)./FErr.^2;

end