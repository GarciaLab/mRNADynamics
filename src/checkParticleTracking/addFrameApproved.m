function [Particles] = addFrameApproved(NChannels, Particles)  
  
  for NCh = 1:NChannels
    if ~isfield(Particles{NCh}, 'FrameApproved')
      for i = 1:length(Particles{NCh})
        Particles{NCh}(i).FrameApproved = true(size(Particles{NCh}(i).Frame));
      end
    else
      for i = 1:length(Particles{NCh})
        if isempty(Particles{NCh}(i).FrameApproved)
          Particles{NCh}(i).FrameApproved = true(size(Particles{NCh}(i).Frame));
        end
      end
    end
  end
  
end