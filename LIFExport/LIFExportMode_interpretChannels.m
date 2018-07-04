function [coatChannel, histoneChannel, fiducialChannel, inputProteinChannel, FrameInfo] = LIFExportMode_interpretChannels(ExperimentType, Channel1, Channel2, Channel3, FrameInfo)
  disp('Interpreting channels...');
  if strcmpi(ExperimentType,'1spot') || strcmpi(ExperimentType,'2spot') || strcmpi(ExperimentType,'2spot1color') || strcmpi(ExperimentType,'inputoutput')
    %Figure out the different channels
    if ~isempty(Channel3)
        Channels={Channel1{1},Channel2{1},Channel3{1}};
    else
        Channels={Channel1{1},Channel2{1}};
    end

    inputProteinChannel = 0; % not used - setting this to 0 to have a common function output among ExperimentTypes

    %Coat protein channel
    coatChannel=find((~cellfun(@isempty,strfind(lower(Channels),'mcp')))|...
        (~cellfun(@isempty,strfind(lower(Channels),'pcp')))|...
        (~cellfun(@isempty,strfind(lower(Channels),'lambda'))));
    if length(coatChannel)>1
        error('Two coat proteins found. Should this be in 2spot2color mode?')
    elseif isempty(coatChannel)    
        error('LIF Mode error: Channel name not recognized. Check MovieDatabase')
    end

    %Histone channel
    histoneChannel=find(~cellfun(@isempty,strfind(lower(Channels),'his')));
    fiducialChannel = histoneChannel;
    %Distinguish between not having histone, but having a dummy channel
    if isempty(fiducialChannel)
        if find(~cellfun(@isempty,strfind(lower(Channels),'dummy')))
            fiducialChannel=0;
        else
            fiducialChannel=0;
            display('Could not find a histone channel. Proceeding without it.')
        end
    end
    % Use MCP-mCherry as a fake histone channel in case there's no
    % Histone channel (Last edited : 3/28/2018, YJK)
    if (fiducialChannel==0)&&...
            ((~isempty(strfind(Channel1{1},'mCherry')))||(~isempty(strfind(Channel2{1},'mCherry'))))
        if (~isempty(strfind(Channel1{1},'mCherry')))
            fiducialChannel=1;
            histoneChannel=1;
        elseif (~isempty(strfind(Channel2{1},'mCherry')))
            fiducialChannel=2;
            histoneChannel=2;
        else
            warning('mCherry channel not found. Cannot generate the fake nuclear image');
        end
    end
      
  elseif strcmpi(ExperimentType,'2spot2color')       %2 spots, 2 colors
    load('ReferenceHist.mat')
    fiducialChannel=0;
    histoneChannel=0;
    coatChannel = 0; % not used - setting this to 0 to have a common function output among ExperimentTypes
    inputProteinChannel = 0; % not used - setting this to 0 to have a common function output among ExperimentTypes

    if (~isempty(strfind(Channel1{1},'mCherry')))||(~isempty(strfind(Channel2{1},'mCherry')))
        if (~isempty(strfind(Channel1{1},'mCherry')))
            fiducialChannel=1;
            histoneChannel=1;
        elseif (~isempty(strfind(Channel2{1},'mCherry')))
            fiducialChannel=2;
            histoneChannel=2;
        else
            warning('mCherry channel not found. Cannot generate the fake nuclear image');
        end
    end
  
  elseif strcmpi(ExperimentType,'input')        %Protein input mode
    %This mode assumes that at least one channel corresponds to the input.
    %It also check whether the second channel is histone. If there is
    %no histone channel it creates a fake channel using one of the
    %inputs.
    
    %Parse the information from the different channels
    Channels = {Channel1{1}, Channel2{1}};
    
    %We have no coat protein here.
    coatChannel = 0;
    
    %Histone channel.
    histoneChannel = find(~cellfun(@isempty,strfind(lower(Channels),'his')));
    if isempty(histoneChannel)
      histoneChannel = 0;
    else
      fiducialChannel = histoneChannel;
    end

    %Input channels
    inputProteinChannel=~cellfun(@isempty,Channels);
    if histoneChannel
        inputProteinChannel(histoneChannel) = 0;
    else
      %If there was no histone channel, we need to choose which
      %input channel to use as our fiducial channel. We'll use
      %the first channel for now. We can try to be smarted about
      %this later on.
      warning('No histone channel found. Finding nuclei using the protein input channel.')
      fiducialChannel = 1;                
    end
    inputProteinChannel = find(inputProteinChannel);
    
    %Save the information about the number of channels in FrameInfo
    for i = 1:length(FrameInfo)
      FrameInfo(i).NChInput = length(inputProteinChannel);
    end
  else
    error('Experiment type not recognized. Check MovieDatabase')
  end
end