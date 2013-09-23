%%                             LoadLsmToMat.m
% Jacques Bothma                            
% Levine Lab, UC Berkeley                       
% Functionally complete                          
%
%
%% Attribution:
%  Feel free to use, modify and distribute this code provided that you
%  attribute Jacques Bothma for development.
%  This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
%  To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/.
%%
% Overview:
%
% Outputs structure containing the slices of z stack with all the image
% info in one slice stored within a single matrix. allows for easier
% manipulation of image data. If there is only one slice the output is a
% simple matrix
%
% Inputs:
%
% filename - string specifying the filename of lsm file contianing the
% image
%
% imnum - stack number to be analyzed.
%
%
%



function  Im = LoadLsmToMat(filename,imnum)

    if exist([filename, '.mat'])~=2                        % Parsing lsm for easy reading
        disp(['Parsing ' filename '.lsm'])
        parselsmwithtime([ filename,'.lsm'])
    else
    end

           
 try       load([filename]);
 catch
     parselsmwithtime([ filename,'.lsm'])
 end
     
      stack = loadlsm(filename,imnum);
        
      
frames= length(stack);


watcher = 0; clear Im;


%%%%%% Store image data in structure

    for f = 1:frames        
        
        SliceBlank=zeros(Datas.LSM_info.DimensionX,Datas.LSM_info.DimensionX,Datas.LSM_info.DimensionChannels, ...
            ['uint' num2str(Datas.Stack1.Image1.TIF.BitsPerSample(1))]);
        
        for i=1:Datas.LSM_info.DimensionChannels
            
            try  % try in case the file does not exist
                
                 SliceBlank(:,:,i) = cell2mat(stack{1,f}(i));
                 
            catch     % stop iterating if exhausted all files.
                disp(['error, cannot find file ', file{i,f}]);
                watcher = 1;
            end
            
            Im.(['Slice' num2str(f)]) = SliceBlank;
        end
          if watcher == 1;
          disp('Loop broken to populate m!!!')
              break; 
          end;      
  
    end
    
    if length(fieldnames(Im))==1
    
        Imt=Im.Slice1;
         clear Im
         Im=Imt;
        
    else
    end
    