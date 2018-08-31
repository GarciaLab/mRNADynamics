function Projection = getHistoneChannelProjection(ProjectionType, HisSlices, ExperimentType, Channel1, Channel2, Channel3, NSlices, ReferenceHist)
  
  if strcmpi(ProjectionType,'medianprojection')
    Projection = median(HisSlices, 3);
  elseif strcmpi(ProjectionType,'middleprojection')
    Projection = max(HisSlices(:,:,round(NSlices*.50):round(NSlices*.75)), [], 3); 
  else
    Projection = max(HisSlices, [], 3);
  end

  %YJK : Think about the case when there is no His channel,
  %and it is inputoutput mode or 1spot mode or 2spot2color.
  %We can use (MCP-mCherry) either inverted or raw
  %images to make fake histone images.
  if (isempty(strfind(Channel1{1}, 'His'))) && (isempty(strfind(Channel2{1}, 'His'))) && (isempty(strfind(Channel3{1}, 'His')))
    if strcmpi(ExperimentType, 'inputoutput')|strcmpi(ExperimentType, '1spot')|strcmpi(ExperimentType,'2spot2color')|strcmpi(ExperimentType,'input')
      if (~isempty(strfind(Channel1{1}, 'NLS')))|(~isempty(strfind(Channel2{1}, 'NLS')))
        %don't invert with NLS-MCP-mCherry
      else
        %We don't want to use all slices. Only the center ones
        StackCenter = round((min(NSlices) - 1) / 2);
        StackRange = StackCenter - 1:StackCenter + 1;
        if strcmp(ProjectionType, 'medianprojection')
            Projection = median(HisSlices(:,:,StackRange), [], 3);
        else
            Projection = max(HisSlices(:,:,StackRange), [], 3);
        end
        %invert images to make nuclei bright
        Projection = imcomplement(Projection);
      end
      Projection = histeq(mat2gray(Projection), ReferenceHist);
      Projection = Projection * 10000;
    end
  end
end
