function idx = find_string_in_cell(myCell,myString)
    % Searches a cell array for a match with a string and returns a logical
    % index for where a match was found
    % Usage:
    % idx = find_string_in_cell(myCell,myString);
    % myCell : a cell array of strings
    % myString : a string
    % 
    % Returns: 
    % idx : logical index where myString matched an entry in myCell
    
    
    idx = find(cellfun(@(x)(...
        ~isempty(strfind(myString,x))),...
        myCell));

end