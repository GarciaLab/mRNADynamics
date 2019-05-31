function [A, B] = addFields(A,B)

%courtesy of Guillaume on MatlabCentral

   for fld = setdiff(fields(A), fields(B))'
      B(end).(fld{1}) = [];   %add missing fields to B
   end
   
   for fld = setdiff(fields(B), fields(A))'
      A(end).(fld{1}) = [];   %add missing fields to A
   end
      
end