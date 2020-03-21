function removeEvenFiles()


for id = 1:length(files)
    
      imFile = '';
      dogStack = imreadStack(imFile);
    
      zSize = dogStack/2;
            
      for zEven = 0:2:zSize
        dogStack(:, :, zEven) = [];
      end
      
      save(imFile, 'dogStack.mat', '-v6');

    
end


end