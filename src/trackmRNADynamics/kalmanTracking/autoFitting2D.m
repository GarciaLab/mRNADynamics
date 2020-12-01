function autoFitting2D()
SpotsIndex = length(Spots{CurrentChannel}(CurrentFrame).Fits)+1;
            breakflag = 0; %this catches when the spot addition was unsuccessful and allows checkparticletracking to keep running and not error out
            use_integral_center = 1;
            
            FitCell = cell(1, liveExperiment.zDim);
            
            for z = 1:liveExperiment.zDim
                spotsIm = imStack(:, :, z);
                try
                    imAbove= imStack(:, :, z+1);
                    imBelow= imStack(:, :, z-1);
                catch
                    imAbove = nan(size(spotsIm,1),size(spotsIm,2));
                    imBelow = nan(size(spotsIm,1),size(spotsIm,2));
                end
                
                Threshold = min(min(spotsIm));
                dog = spotsIm;
                im_thresh = dog >= Threshold;
                [im_label, ~] = bwlabel(im_thresh);
                microscope = FrameInfo(1).FileMode;
                show_status = 0;
                fig = [];
                k = 1; %This is supposed to be the index for the particles in an image.
                %However, this image only contains one particle
                neighborhood_px = round(1300 / liveExperiment.pixelSize_nm); %nm
                %Get the information about the spot on this z-slice                
                Fit = identifySingleSpot(k, {spotsIm,imAbove,imBelow}, im_label, dog, neighborhood_px, liveExperiment.snippetSize_px, ...
                    liveExperiment.pixelSize_nm, show_status, fig, microscope,...
                    [1, ConnectPositionx, ConnectPositiony], [], '', CurrentFrame, [], z);
             
                
                FitCell{z} = Fit;
                
                
            end
            Fits = [];
            
            for z = 1:liveExperiment.zDim
                if ~isempty(FitCell{z})
                    fieldnames = fields(FitCell{z});
                    if isempty(Fits)
                        Fits = FitCell{z};
                    else
                        for field = 1:length(fieldnames)
                            Fits.(fieldnames{field}) =  [Fits.(fieldnames{field}), FitCell{z}.(fieldnames{field})];
                        end
                    end
                end
            end
            
            force_z = 0;
            [Fits,~, ~] = findBrightestZ(Fits, -1, use_integral_center, force_z, []);