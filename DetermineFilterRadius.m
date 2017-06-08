function OptimalRadius = DetermineFilterRadius(NucEven,StartRadius,EndRadius,NumberRadiusPoints,scalingfactor,numberhistogrampoints,Offset)
% OptimalRadius = DetermineFilterRadius(NucEven,StartRadius,EndRadius,
%           NumberRadiusPoints,scalingfactor,numberhistogrampoints,Offset)
%
% DESCRIPTION
% This code determines the optimal size for the difference of gaussian
% filter that does the nuclear segmentation. This is important when the
% size of the nuclei vary i.e one is analysing embryos at a number of
% different cell cycles. It looks at how the number of objects in an image
% changes with the radius of the LOG filter. At the
% point where the radius approaches the optimal value the number of objects
% as a function of radius will plateau. The middle point of this plateau 
% plus some offset is used as the optimal radius.
%
% PARAMETERS
% NucEven: Evening filtered image of the nuclei
% StartRadius: Smallest radius of gaussian filter to smooth nuclear image 
%              (~40)
% EndRadius: Largest radius of gaussian filter to smooth nuclear image
%            (~20)
% NumberRadiusPoints: Number of radius vlaues between StartRadius and
%                     EndRadius to use.
% scalingfactor: Scalingfactor to use reduce size of image to speed up
%                processing.
% numberhistogrampoints: Related to how wide a region of the pateau to use
%                        for determining optimal value.
% Offset: Extra size to add to radius (empirical)
%
% OPTIONS
% None
%
% CONTROLS
% None
%
% OUTPUT
% OptimalRadius: Optimal filter radius
%
% Author (contact): Unknown
% Created: Unknown
% Last Updated: Unknown
%
% Commented by: Meghan Turner (meghan_turner@berkeley.edu), 6/8/17

troubleshooting = 0;

[ro,co]=size(NucEven);     % Resizing of the image to speed up calculations
r = round(scalingfactor*ro);
c = round(scalingfactor*co);
NucS=imresize(NucEven,[r,c]);    

RadiusRange = linspace(StartRadius,EndRadius,NumberRadiusPoints);  % Defin-
                                                 % ing the different filter
                                                 % radius values to use

NumObj=zeros(NumberRadiusPoints,2);         % Initializing array that will 
                                            % store the number of objects 
                                            % at a different radius

for i=1:NumberRadiusPoints           % Determining the number of objects at
                                     % the different radius values
bw = segmentnucleiCoreNoWatershed(NucS,round(10*RadiusRange(i)*scalingfactor),RadiusRange(i)*scalingfactor,0);
NumObj(i,1) = max(max(bwlabeln(bw)));
NumObj(i,2) =RadiusRange(i);
end

if troubleshooting                   % Plots the number of objects as a 
                                     % function of filter radius during 
                                     % troubleshooting
figure, plot(NumObj(:,2),NumObj(:,1),'.r')
title('Number of objects as a function of filter radius')
end

[FreqNumObj,NumObjInt] = hist(NumObj(:,1),numberhistogrampoints);  % Deter-
                                            % mines the histogram of the 
                                            % number of objects used to 
                                            % locate plateau

if troubleshooting
figure, bar(NumObjInt,FreqNumObj)
title('Histogram of the number of objects to find plateau')
end

NumObjIntatMaxFreq = NumObjInt(find(FreqNumObj==max(FreqNumObj))); %Deter-
                                            % mines the most often occuring
                                            % number of objects

if length(NumObjIntatMaxFreq)>1;
    NumObjIntatMaxFreq=mean(NumObjIntatMaxFreq);
end

HalfWidth = abs(NumObjInt(1)-NumObjInt(2))/2 ; %HalfWidth of the bins


PlateauValues = RadiusRange(NumObj(:,1)>=NumObjIntatMaxFreq-HalfWidth & NumObj(:,1)<=NumObjIntatMaxFreq+HalfWidth); %The radius values that are on the plateau

OptimalRadius=median(PlateauValues)+Offset; %Determining the optimal Radius


