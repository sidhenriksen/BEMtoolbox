function matFiles = get_mat_files(dataDir)
    
    allFiles = dir(dataDir);
    allFileNames = {allFiles.name};
    mat_strcmp = @(x)(any(strfind(x,'.mat')));
    matFilesIdx = cellfun(mat_strcmp,allFileNames);
    
    matFiles = allFileNames(matFilesIdx);
end
