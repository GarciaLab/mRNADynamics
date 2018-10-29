function Projection = getDefaultChannelProjection(ProjectionType, NSlices, HisSlices)
  %We don't want to use all slices. Only the center ones
  StackCenter = round((min(NSlices) - 1) / 2);
  StackRange = StackCenter - 1:StackCenter + 1;
  otherSlices = HisSlices;
  if strcmp(ProjectionType, 'medianprojection')
    Projection = median(HisSlices(:,:,StackRange), 3);
    otherProjection = median(otherSlices(:,:,StackRange), 3);
  else
    Projection = max(HisSlices(:,:,StackRange), [], 3);
    otherProjection = max(otherSlices(:,:,StackRange), [], 3);
  end
  Projection = Projection + otherProjection;
end
