function parentFolder = parentFolder(folder)

mydir  = folder;
idcs   = [strfind(mydir,'/'), strfind(mydir,'\')];
parentFolder = mydir(1:idcs(end)-1);

end