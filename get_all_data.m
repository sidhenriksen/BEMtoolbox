function [allCell,allModel] = get_all_data(runData)
    % Runs or loads the data from different analyses, 
    % merges and returns the corresponding data structs.
    % Usage:
    % [allCell, allModel] = get_all_data(runData)
    % runData : 
    %               0 - load from disk, 
    %               1 - run sims, overwrite 
    %               2 - run sims, do not overwrite     

    saveData = 1; 
    
    if isunix
        dataDir = '/sid/nim/data/';
        outDataDir = '/sid/nim/analysis/fit_data/';
    else
        dataDir = 'X:/readlab/Sid/nim/data/';
        outDataDir = 'X:/readlab/Sid/nim/analysis/fit_data/';
    end
    
    
    
    matFiles = get_mat_files(dataDir);        
    
    contrastData = struct();
    
    twopassData = struct();

    for k = 1:length(matFiles);

        fileName = matFiles{k};

        outFileName = [outDataDir,fileName(1:end-4),'_fit.mat'];                                        

        if runData % if this is true, we want to actually run the computations

            data=load(fileName);

            % if the cross-validation has been run, use the determined best number
            % of subunits; otherwise use default settings
            if isfield(data,'nExc');

                nExc = data.nExc;

                nInh = data.nInh;

            else

                nExc = 3;

                nInh = 1;

            end
            
            % Skip if outFileName already exists and runData says to do
            % so.
            if exist(outFileName,'file') && runData == 2; continue; end
                        

            NimCell = data.NimCell;     
            

            NimModel = fit_NimModel(NimCell,nExc,nInh);
            
            NimModel.rEst = poissrnd(NimModel.rEst);
            
            NimCell.NimFit = NimModel.NimFit; % get_contrast_data needs this to run
            
            contrastCell = get_contrast_data(NimCell);
            
            contrastModel = get_contrast_data(NimModel);
            
            
            twopassCell = get_twopass_data(NimCell);
            
            twopassModel = get_twopass_data(NimModel);            
            
                        
            if saveData

                save(outFileName,'twopassCell','twopassModel','contrastCell','contrastModel');

            end

        else % if runData=0, then simply load the precomputed data for this cell

            load(outFileName);

        end
        
        
        contrastData(k).Cell = contrastCell;

        contrastData(k).Model = contrastModel;
                

        twopassData(k).Cell = twopassCell;        

        twopassData(k).Model = twopassModel;
        
        twopassData(k).Cell.fileName = fileName;
        twopassData(k).Model.fileName = fileName;

    end
        
    
        
    
    
    [allCell,allModel] = custom_merge_structs(twopassData,contrastData);
    

    
    
end

function [allCell,allModel] = custom_merge_structs(twopassData,contrastData);

    assert(length(twopassData)==length(contrastData),'Must be same length');
    
    %allCell = struct();
    %allModel = struct();
        
    for k = 1:length(twopassData);
        cellData = merge_structs(twopassData(k).Cell,contrastData(k).Cell);
        modelData = merge_structs(twopassData(k).Model,contrastData(k).Model);
        
        
            
        allCell(k) = cellData;
        allModel(k) = modelData;
    end

end
