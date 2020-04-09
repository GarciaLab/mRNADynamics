function myCleanupFun()
warning('off',   'MATLAB:class:DestructorError');
delete(findall(0));
    
end