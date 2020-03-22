function addJavaPathsForLivemRNA()

if ~exist([userpath, filesep, 'javaclasspath.txt'], 'file')
    path1 = 'C:\Program Files\Weka-3-8-4\weka.jar';
    path2 = 'C:\Users\Armando\Desktop\fast random forest\fastrandomforest-2019.12.3.jar';
    path3 = 'C:\Users\Armando\Downloads\Fiji.app\jars\imagescience-3.0.0.jar';

    fileID = fopen([userpath, filesep, 'javaclasspath.txt'],'w');

    fprintf(fileID,'%s \n %s \n %s', path1, path2, path3);
    
    fclose(fileID);

end
