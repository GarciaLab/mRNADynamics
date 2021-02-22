% Vish:
% cast a variable to specified type. 
% Will cast to guint8 if variable is GPU-type, and
% uint8 if variable is CPU type.
function a = smartcast(a,prec)

  isjacket = 1;
  if isempty(regexp(path,'[jJ]acket[/\\]engine'))
    isjacket = 0;
  end

  ty = isa(a,'garray');
  if ty
    a = feval(['g' prec], a);
    return
  else
    a = feval([prec], a);
    return
  end
end
