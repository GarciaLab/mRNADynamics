function keyInputHandler = NuclearTrackingEventHandler(cptState)
 
    function keyInput(cc)
        if cc == 'l'
            % Split a nucleus and select one or two daughter cells or stop the lineage
            warning('NL: code does not currently support splitting schnitz traces.')
%             [cptState.Particles, cptState.PreviousParticle, cptState.schnitzcells] = splitNuclei(cptState.schnitzcells, ...
%                 cptState.CurrentFrame, cptState.CurrentChannel, cptState.CurrentParticle, cptState.Particles);
        elseif cc == '2'
            % Set parent of current nucleus
            cptState.schnitzcells = setParentNucleus(cptState.schnitzcells, ...
                cptState.CurrentFrame, cptState.CurrentChannel, cptState.CurrentParticle, cptState.Particles);
        elseif cc == '$'
            % Add particle to nucleus
            cptState.Particles = addNucleusToParticle(cptState.Particles, cptState.CurrentFrame, ...
                cptState.CurrentChannelIndex, cptState.UseHistoneOverlay,...
                cptState.schnitzcells, cptState.CurrentParticle);
        end
    end

    keyInputHandler = @keyInput;
end
