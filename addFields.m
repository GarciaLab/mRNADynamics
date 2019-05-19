function [A, B] = addFields(A,B)

%courtesy of Guillaume on MatlabCentral

   for fld = setdiff(fieldnames(A), fieldnames(B))'
      A(end).(fld{1}) = [];   %add missing fields to A
   end
   
   for fld = setdiff(fieldnames(B), fieldnames(A))'
      B(end).(fld{1}) = [];   %add missing fields to B
   end
      
end