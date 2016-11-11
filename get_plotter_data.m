function [allCell,allModel] = get_plotter_data()
    % Loads the data needed for the plotter
        
    
    dataDir = get_data_dir();
    
    matFiles = get_mat_files(dataDir);        
    
    contrastData = struct();
    
    twopassData = struct();

    for k = 1:length(matFiles);

        load([dataDir,'/',matFiles{k}]);

        contrastData(k).Cell = contrastCell;

        contrastData(k).Model = contrastModel;
                

        twopassData(k).Cell = twopassCell;        

        twopassData(k).Model = twopassModel;
        
        fileName = get_file_name(matFiles{k});
        
        twopassData(k).Cell.fileName = fileName;
        twopassData(k).Model.fileName = fileName;

    end
        
    
        
    
    
    [allCell,allModel] = custom_merge_structs(twopassData,contrastData);
    

    
    
end

function fileName = get_file_name(matFile);

    underscore_idx = strfind(matFile,'_');
    
    fileName = [matFile(1:(underscore_idx-1)),'.mat'];

end

function [allCell,allModel] = custom_merge_structs(twopassData,contrastData);

    assert(length(twopassData)==length(contrastData),'Must be same length');
        
        
    for k = 1:length(twopassData);
        cellData = merge_structs(twopassData(k).Cell,contrastData(k).Cell);
        modelData = merge_structs(twopassData(k).Model,contrastData(k).Model);                
            
        allCell(k) = cellData;
        allModel(k) = modelData;
    end

end


function dataDir = get_data_dir();
    
    fullpath = mfilename('fullpath');
    
    if isunix
        slash = '/';
    else
        slash = '\';
    end
    
    slashIndices = findstr(fullpath,slash);
    
    lastIdx = slashIndices(end);
    
    dataDir = [fullpath(1:lastIdx),'fit_data'];
end