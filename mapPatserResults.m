function mapPatserResults(path)
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

    close all force;
    t = fileread(path);
    tinverse = flip(t);
    cutOffExpression = '(Print scores greater than or equal to \=?)(\d*\.\d)';
    cutOff = regexp(t,cutOffExpression,'tokens');
    cutOffScore = str2num(char(cutOff{1}(2)));
    
    
    expression0 = '(\=erocs)';
    expression1 = 'position='; % to retrieve the location information
    s0 = regexp(tinverse, expression0);
    s1 = regexp(t, expression1); % locating all position values
    tf = {};
    s02 = cell(length(s0),1);
    for i = 1:length(s0)
        pos0 = s0(i);
        postf = strfind(tinverse(pos0:length(tinverse)), '_FT');
        
        s02{i} = t((s1(i)+length(expression1)):(s1(i)+length(expression1)+8)); % 8 = # of characters before score=
%         try
%             postf = strfind(tinverse(pos:length(tinverse)), '_FT');
%         catch
%             postf = strfind(tinverse(pos:length(tinverse)), '_FT');  
%         end
        b= postf(1)+pos0;
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
    siteIndices = [];
    complementLocation = [];
    complementSiteIndices = [];
    for i = 1:length(s01)
        s = s01{i};
        s2 = s02(i); 
        values(i) = str2double(flip(s{1})); % flip does not do anything?
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
    try
        standardizeFigure(ax, [], 'bar', b, 'red');
    end

    
    figure(2)%,'units','normalized','outerposition',[0 0 1 1]);
    % Sites versus Position
    clf
    hold on
    
    displayRangeForward = linspace(cutOffScore,max(values),length(location));  %transpose([0;0;1]*linspace(0.4,1,length(location)));   
    scatter(location,-1*ones(1,length(location)),[],displayRangeForward,'filled'); % locations
    
    displayRangeReverse = linspace(cutOffScore,max(values),length(complementLocation));%transpose([0;0;1]*linspace(0.4,1,length(complementLocation)));    
    scatter(complementLocation,ones(1,length(complementLocation)),[],displayRangeReverse,'filled');
    
    ax = gca;
    ylim(ax,[-2,2]);
    ax.YTick = -2:2;
    set(ax,'YTickLabel',{' ','forward',' ', 'reverse',' '});
    ylabel('orientation')
    xlabel('location')
    backslashes = strfind(path,'\');
    title({path(backslashes(length(backslashes))+1:length(path)-4), ['cut off score: ' num2str(cutOffScore)]},'Interpreter','none')
    h = colorbar;
    h.Limits = [cutOffScore ceil(max(values))];
    
    % sorting by location for labeling 
    [location,siteIndices] = sortLocations(location,siteIndices);
    [complementLocation,complementSiteIndices] = sortLocations(complementLocation,complementSiteIndices);
     
    arrowAdjustment = 0;
    for i = 1:length(location)
        origin = -1;
        arrowPosition = [origin*1.15,origin*0.85];
        wordPosition = [origin*1.25,origin*0.75];
        rotation = [-45,45];
        arrow = {'$\uparrow$','$\downarrow$'};
        text(location(i)-arrowAdjustment,arrowPosition(1+mod(i,2)),arrow{1+mod(i,2)},'Interpreter','latex');
        text(location(i),wordPosition(1+mod(i,2)),strcat(tf(siteIndices(i))),'Rotation',rotation(1+mod(i,2)),'Interpreter','none');
    end 
    
    for i = 1:length(complementLocation)
        origin = 1;
        arrowPosition = [origin*1.15,origin*0.85];
        wordPosition = [origin*1.25,origin*0.75];
        rotation = [45,-45];
        arrow = {'$\downarrow$','$\uparrow$'};
        
        text(complementLocation(i)-arrowAdjustment,arrowPosition(1+mod(i,2)),arrow{1+mod(i,2)},'Interpreter','latex');
        text(complementLocation(i),wordPosition(1+mod(i,2)),strcat(tf(complementSiteIndices(i))) ,'Rotation',rotation(1+mod(i,2)),'Interpreter','none');
    end
    
    try
        standardizeFigure(ax, [], 'axeslinewidth', 1);
    end
    ax.XGrid = 'on';
    ax.XMinorGrid = 'on';
    
    hold off

%% Saving Compiled Information    
%     fileName = [path(1:backslashes(length(backslashes))) 'condensedPatserResult_' path(backslashes(length(backslashes))+1:length(path))];
%     fileID = fopen(fileName,'w');
%     for i = 1:length(values)
%         fprintf(fileID,'%6s %12s\n','x','exp(x)');
%     end
%     
%     fclose(fileID);

    
    
    
end

function [sortedLocations,sortedIndices] = sortLocations(locations,indices)
    
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
