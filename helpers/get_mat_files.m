function matFiles = get_mat_files(dataDir,skipEnergyModel)
    
    if nargin < 2
        skipEnergyModel = false;
    end
    
    allFiles = dir(dataDir);
    allFileNames = {allFiles.name};
    mat_strcmp = @(x)(any(strfind(x,'.mat')));
    matFilesIdx = cellfun(mat_strcmp,allFileNames);
    
    matFiles = allFileNames(matFilesIdx);
    
    if skipEnergyModel
        
        matFiles = matFiles(~find_string_in_cell(matFiles,'energymodel.mat',1));
        
    end
end
