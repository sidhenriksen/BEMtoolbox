function modelName = get_model_file(cellFile)
    
    patternToMatch = '[a-z][a-z][a-z][A-Z][0-9][0-9][0-9]c[0-9]*.mat';
    
    idx = regexp(cellFile,patternToMatch);
    
    modelName = cellFile(idx:end);
    
    modelName(idx+7) = 'm';
    

end

