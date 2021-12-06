function [Ellipses, schnitzcells] = addStrayEllipsesToSchnitzcells(Ellipses, schnitzcells)




schnitzcellsOld = schnitzcells;


for f = 1:length(Ellipses)
   
   error_message = 'Ellipses structure missing a column. Did you accidentally select ''noTrack''?';
   assert( size(Ellipses{f}, 2) > 8, error_message)

    
   strayEllipseIndexes = find(Ellipses{f}(:, 9) == 0);
   
   for s = 1:strayEllipseIndexes
       %create new schnitz and add ellipse reference
       schnitzcells(end+1).cellno = s;
       
       %fill up the other fields
       schnitzcells(end).frames = uint16(f);
       schnitzcells(end).cenx = Ellipses{f}(s, 1);
       schnitzcells(end).ceny = Ellipses{f}(s, 2);
       schnitzcells(end).len = single(mean([Ellipses{f}(s, 3), Ellipses{f}(s, 4)]));
       
       
   end

      
end


%some validation
assert(length(schnitzcellsOld) <= length(schnitzcells));

end