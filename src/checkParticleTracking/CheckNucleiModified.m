function CheckNucleiModified(cptState, DropboxFolder, Prefix, fish)

	if cptState.nucleiModified
	    Ellipses = cptState.Ellipses; % assign it to a local variable so we can save it to its .mat below
	    save([DropboxFolder, filesep, Prefix, filesep, 'Ellipses.mat'], 'Ellipses');
	    
	    % Decide whether we need to re-track
	    userPrompt = 'Did you make changes to nuclei and thus require re-tracking? (y/n)';
	    reTrackAnswer = inputdlg(userPrompt);
	    if contains(reTrackAnswer,'n')
	        disp('Ellipses saved. Per user input, not re-tracking. Exiting.')
	    else
	        opts = {};
	        if fish
	            opts = [opts, 'markandfind'];
	        end
	        disp('Ellipses saved. Running TrackNuclei to incorporate changes.')
	        TrackNuclei(Prefix,'NoBulkShift','ExpandedSpaceTolerance', 1.5, 'retrack', 'nWorkers', 1, opts{:});
	    end
	end

end
