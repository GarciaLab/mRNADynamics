function Projections = calculateOMETIFFProjections(ProjectionType, NuclearPlanes)
  % load the reference histogram used to enhance contrast on histone channel
  load('ReferenceHist.mat');
  
  
  
  % NuclearPlanes(CurrentZSlize, CurrentNFrame) = TIFImages{1}{i};
  
  NFrames = size(NuclearPlanes);
  NFrames = NFrames(2);
  
  for i = 1:NFrames
    
    % pre-allocate projections for current frame
    zSize = size(NuclearPlanes,1);
    HisSlices = uint16(zeros([size(NuclearPlanes{1}, 1), size(NuclearPlanes{1}, 2), zSize]));
    
    for z = 1:zSize
      zSlice = NuclearPlanes(:,i);
      zSlice = zSlice(z);
      HisSlices(:, :, z) = zSlice{1};
    end
  
    if strcmpi(ProjectionType, 'medianprojection')
      Projection = median(HisSlices, 3);
    elseif strcmpi(ProjectionType, 'middleprojection')
      Projection = max(HisSlices(:, :, round(zSize * .50):round(zSize * .75)), [], 3);
    elseif strcmpi(ProjectionType, 'maxprojection')
      Projection = max(HisSlices, [], 3);
    else
      ProjectionBounds = strsplit(ProjectionType, ':');
      SortedHisSlices = sort(HisSlices, 3, 'descend');
      max_custom = str2double(ProjectionBounds{2});
      min_custom = str2double(ProjectionBounds{3});
      Projection = mean(SortedHisSlices(:, :, max_custom:min_custom), 3);
    end

    Projection = histeq(mat2gray(Projection), ReferenceHist);
    Projection = Projection * 10000;
    
    Projections{i} = Projection;
  end
  
end