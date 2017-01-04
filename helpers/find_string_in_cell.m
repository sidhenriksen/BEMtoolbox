function idx = find_string_in_cell(myCell,myString,returnLogical)
    % Searches a cell array for a match with a string and returns a logical
    % index for where a match was found
    % Usage:
    % idx = find_string_in_cell(myCell,myString);
    % myCell : a cell array of strings
    % myString : a string
    % returnLogical : boolean whether to return a logical or index array
    % (default is index)
    % 
    % Returns: 
    % idx : logical index where myString matched an entry in myCell
    
    if nargin < 3
        returnLogical = 0;
    end
    
    
    idx = cellfun(@(x)(...
        ~isempty(strfind(myString,x))),...
        myCell);
    
    if ~returnLogical
        idx = find(idx);
    end

end