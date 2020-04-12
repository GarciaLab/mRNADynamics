function mapPatserResults(path, varargin)
% mapPatserResults(path)
%
% DESCRIPTION
% Maps Patser scores onto the DNA sequence of interest.
%
% ARGUMENTS
% path: Path to your Paster result text file. Save the output from http://stormo.wustl.edu/consensus/cgi-bin/Server/Interface/patser.cgi
%
% OPTIONS
% None.

% OUTPUT
% A bar chart showing relevant binding sites and their scores, and a map of
% the locus showing relevant binding site locations along with their
% scores.
%
% Author (contact): Armando Reimer & Emma Luu (areimer@berkeley.edu)
% Created:
% Last Updated: 11/29/2017
%
% Documented by: Emma Luu(emma_luu@berkeley.edu)

offset = [];

for arg = 1:length(varargin)
    if strcmpi(varargin{arg}, 'index')
        offset = varargin{arg+1};
    end
end

close all force;
t = fileread(path);
tinverse = flip(t); %flipping text for score search
cutOffExpression = '(Print scores greater than or equal to \=?)(\d*\.\d)';
cutOff = regexp(t,cutOffExpression,'tokens');
cutOffScore = str2double(char(cutOff{1}(2)));


expression0 = '(\=erocs)';
expression1 = 'position='; % to retrieve the location information
s0 = regexp(tinverse, expression0);
s1 = regexp(t, expression1); % locating all position values
tf = {};
s02 = cell(length(s0),1);
for i = 1:length(s0)
    pos0 = s0(i);
    postf = strfind(tinverse(pos0:length(tinverse)), '_FT');
    
    poswidth = strfind(tinverse(pos0:length(tinverse)), ' :x');
    s02{i} = t((s1(i)+length(expression1)):(s1(i)+length(expression1)+8)); % 8 = # of characters before score=
    b= postf(1)+pos0;
    c = poswidth(1)+pos0;
    width{i} = flip(tinverse(c-3:c-2));
    tf{i} = tinverse(b-50:b-2);
    z = strfind(tf{i}, sprintf('\r'));
    tf{i} = tf{i}(z+4:end);
    tf{i} = flip(tf{i});
end
expression01 = '(\d*\.\d)(?=   \=erocs)';
s01 = regexp(tinverse, expression01, 'tokens');
values = zeros(1,length(s01));
complementPositionCounter = 0; % on the reverse orientation
positionCounter = 0; %on the forward orentiation
location = [];
footprints = cellfun(@str2double, width);
siteIndices = [];
complementLocation = [];
complementSiteIndices = [];
for i = 1:length(s01)
    s = s01{i};
    s2 = s02(i);
    values(i) = str2double(flip(s{1}));
    % making different arrays for sequence and complement sequence locations.
    if isnan(str2double(s2)) %checking if this the complement (#C)
        complementPositionCounter = complementPositionCounter + 1;
        complementLocation(complementPositionCounter) = str2double(erase(s2,'C')); % Removing the C
        complementSiteIndices(complementPositionCounter) = i;% an array to store the indices of the sites.
    else
        positionCounter = positionCounter + 1;
        location(positionCounter) = str2double(s2);
        siteIndices(positionCounter) = i;% an array to store the indices of the sites.
    end
end
tf = flip(tf);
values = flip(values);
footprints = flip(footprints);


%sorting by score for colorbar
[location,siteIndices] = sortScore(location,values,siteIndices);
[complementLocation,complementSiteIndices] = sortScore(complementLocation,values,complementSiteIndices);


%% Plotting
figure(1) % Scores versus Sites
b = bar(values);
ax = gca;
set(ax,'xticklabel',tf)
ax.XTick = 1:1:length(tf);
ax.XTickLabelRotation=45;
ax.TickLabelInterpreter = 'none';
xlabel('sites')
ylabel('scores')
backslashes = strfind(path,'\');
title(path(backslashes(length(backslashes))+1:length(path)-4),'Interpreter','none')
standardizeFigure(ax, [], 'red');


figure(2); % Sites versus Position
clf
hold on
% sorting by location for graph
[location,siteIndices] = sortLocations(location,siteIndices);
[complementLocation,complementSiteIndices] = sortLocations(complementLocation,complementSiteIndices);

orientationAxisAdjustment = [-1 1]*(1/8);
orientationAxisShift = repmat(orientationAxisAdjustment,1,...
    max(length(location),length(complementLocation)));%Alternating for increased visibility
markerSize = 36*5; % Default = 36
separationY = 2; % This will change the distance of the orientations

scatter(location,-separationY*ones(1,length(location))+ orientationAxisShift(1:length(location)),...
    markerSize,values(siteIndices),'filled');%displayRangeForward,'filled'); % locations
scatter(complementLocation,separationY*ones(1,length(complementLocation))+ orientationAxisShift(1:length(complementLocation)),...
    markerSize,values(complementSiteIndices),'filled');

ax = gca;
ylim(ax,[-2 2]*separationY);
ax.YTick = [-2:2]*separationY;
set(ax,'YTickLabel',{' ','Forward',' ', 'Reverse',' '});
ax.XGrid = 'on';
ax.XMinorGrid = 'on';
ylabel('Orientation','FontSize',20,'FontWeight','bold')
xlabel('Location','FontSize',20,'FontWeight','bold')
backslashes = strfind(path,'\');
title({path(backslashes(length(backslashes))+1:length(path)-4), ['Cut Off Score: ' num2str(cutOffScore)]},'Interpreter','none','FontSize',24,'FontWeight','bold')
h = colorbar;
h.Limits = [cutOffScore ceil(max(values))];
h.Label.String = 'Score (AU)';
h.Label.FontSize = 20;

arrowXAdjustment = -0.5; % This is for a horizontal shift for the arrow
arrowPercentageShift = 0.15; % This is to shift the arrow's position
wordPercentageShift = 0.25; % This is to shift the word's position
%adding binding site labels for the forward orientation
for i = 1:length(location)
    origin = -1*separationY;
    upShift = 1 + [arrowPercentageShift, wordPercentageShift];
    downShift = 1 - [arrowPercentageShift, wordPercentageShift];
    arrowPosition = [origin*downShift(1),origin*upShift(1)]+ flip(orientationAxisAdjustment);
    wordPosition = [origin*downShift(2),origin*upShift(2)]+ flip(orientationAxisAdjustment);
    rotation = [45,-45];
    arrow = {'$\downarrow$','$\uparrow$'};
    text(location(i)+arrowXAdjustment,arrowPosition(1+mod(i,2)),arrow{1+mod(i,2)},'Interpreter','latex');
    text(location(i),wordPosition(1+mod(i,2)),strcat(tf(siteIndices(i))),'Rotation',rotation(1+mod(i,2)),'Interpreter','none', 'FontSize', 15, 'FontWeight', 'bold');
end

%adding binding site labels for the reverse orientation
for i = 1:length(complementLocation)
    origin = separationY;
    arrowPosition = [origin*1.15,origin*0.85]+flip(orientationAxisAdjustment);
    wordPosition = [origin*1.25,origin*0.75]+flip(orientationAxisAdjustment);
    rotation = [45,-45];
    arrow = {'$\downarrow$','$\uparrow$'};
    
    text(complementLocation(i)+arrowXAdjustment,arrowPosition(1+mod(i,2)),arrow{1+mod(i,2)},'Interpreter','latex');
    text(complementLocation(i),wordPosition(1+mod(i,2)),strcat(tf(complementSiteIndices(i))) ,'Rotation',rotation(1+mod(i,2)),'Interpreter','none', 'FontSize', 15, 'FontWeight', 'bold');
end

%standardizeFigure(ax, [], 'axeslinewidth', 1);
ax.XGrid = 'on';
ax.XMinorGrid = 'on';

hold off

%% Saving Compiled Information
textFileName = 'CondensedResults.txt';

%Saving the condensed results in the same location as the given path.
fileName = [path(1:length(path)-4) textFileName];
fileID = fopen(fileName,'w');
fprintf(fileID, 'Forward Direction:\r\n');
fprintf(fileID, '%16s %9s %10s \r\n','Transcription F.', 'Score', 'Location');
for i = 1:length(location)
    fprintf(fileID, '%16s %9g %10g\r\n', char(tf(siteIndices(i))),...
        values(siteIndices(i)),location(i));
end

fprintf(fileID, '\r\nReverse Direction:\r\n');
fprintf(fileID, '%16s %9s %10s \r\n','Transcription F.', 'Score', 'Location');
for i = 1:length(complementLocation)
    fprintf(fileID, '%16s %9g %10g\r\n', char(tf(complementSiteIndices(i))),...
        values(complementSiteIndices(i)),complementLocation(i));
end
fclose(fileID);

if ~isempty(offset)
%%
    %make a .gb file with annotations

    fileName = [path(1:length(path)-4),'.gb'];
    fid = fopen(fileName, 'w');
    fidIn = fopen('E:\Armando\LivemRNA\Data\Dropbox\pib-phsp70-ms2v5-lacz.gb', 'r');


    fin = fscanf(fidIn, '%c');
    qI = strfind(fin, 'Qualifiers');
    fprintf(fid, '%c', fin(1:qI+length('Qualifiers')));


    offsetLocation = location + offset - 1;
    offsetComplementLocation = complementLocation + offset - 1;

    for i = 1:length(offsetLocation)

        endLocation = offsetLocation(i) + footprints(siteIndices(i)) - 1;

        if i == 1
            fprintf(fid, '%c\r\n', '');
        end
        fprintf(fid, '%c', ['     misc_feature    ',num2str(offsetLocation(i)), '..', num2str(endLocation)]);
        fprintf(fid, '%c\r\n', '');
        fprintf(fid, '%c', ['                     /label="',tf{siteIndices(i)},' score: ',num2str(round(values(siteIndices(i)))),'"']);
        if i ~= length(offsetLocation)
            fprintf(fid, '%c\r\n', '');
        end
    end


    %not sure if the complement location is at the end or beginning
    for i = 1:length(offsetComplementLocation)

        endLocation = offsetComplementLocation(i) + footprints(complementSiteIndices(i)) - 1;

        if i == 1
            fprintf(fid, '%c\r\n', '');
        end
        fprintf(fid, '%c', ['     misc_feature    complement(',num2str(offsetComplementLocation(i)), '..', num2str(endLocation),')']);
        fprintf(fid, '%c\r\n', '');
        fprintf(fid, '%c', ['                     /label="',tf{complementSiteIndices(i)},' score: ',num2str(round(values(complementSiteIndices(i)))),'"']);
        fprintf(fid, '%c\r\n', '');
    end
    fprintf(fid, '%c', fin(qI+length('Qualifiers')+1:end));
    fclose(fid);
end

end

function [sortedLocations,sortedIndices] = sortLocations(locations,indices)
% This function will sort the sites by locations (in increasing order).
for i = 1:(length(locations)-1)
    min = locations(i);
    minIndex = i;
    
    for j = i+1 : length(locations)
        if locations(j) < min
            min = locations(j);
            minIndex = j;
        end
    end
    locations(minIndex) = locations(i);
    locations(i) = min;
    temp = indices(minIndex);
    indices(minIndex) = indices(i);
    indices(i) = temp;
end
sortedLocations = locations;
sortedIndices = indices;
end


function [sortedLocations,sortedIndices] = sortScore(locations, values,indices)
% This function sorts the sites by scores (in increasing order).
for i = 1:(length(locations)-1)
    minValue = values(indices(i));
    minLocation = locations(i); % the corresponding location to the minValue
    minIndex = i;
    
    for j = i+1 : length(locations)
        if values(indices(j)) < minValue
            minValue = values(indices(j));
            minLocation = locations(j);
            minIndex = j;
        end
    end
    locations(minIndex) = locations(i);
    locations(i) = minLocation;
    temp = indices(minIndex);
    indices(minIndex) = indices(i);
    indices(i) = temp;
end
sortedLocations = locations;
sortedIndices = indices;
end
