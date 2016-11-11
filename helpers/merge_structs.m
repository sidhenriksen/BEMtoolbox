function mergedStruct = merge_structs(struct1,struct2,varargin)
    % Function to merge two or more structures.
    % Usage:
    % mergeStruct = merge_structs(struct1,struct2,...)
    %
    % Note: the structures must have equal length
        
    assert(length(struct1)==length(struct2),'AssertionError: Lengths of the two structs must be the same');
        
    mergedStruct = struct_merge(struct1,struct2);    
    
    % call this recursively to merge all structs
    if nargin == 3
        
        mergedStruct = merge_structs(mergedStruct,varargin{1});
        
    elseif nargin > 3
        
        mergedStruct = merge_structs(mergedStruct,varargin{1},varargin(2:end));
        
    end
    
end

function myStruct = struct_merge(struct1,struct2)

    myStruct = struct();
    
    struct1Fields = fields(struct1);
    
    struct2Fields = fields(struct2);
        
    mergeStruct1Fields = struct1Fields;
    
    mergeStruct2Fields = struct2Fields;
    
    intersectFields = intersect(struct1Fields,struct2Fields);
    
    if ~isempty(intersectFields)
        
        for k = 1:length(intersectFields);
           
            
           
            idx1 = find(strcmp(struct1Fields,intersectFields{k}));
            
            idx2 = find(strcmp(struct2Fields,intersectFields{k}));
            
            v1 = struct1.(intersectFields{k});
            v2 = struct2.(intersectFields{k});
            
            if ~issame(v1,v2);
                      
                mergeStruct1Fields{idx1} = [intersectFields{idx1},'_s1'];
            
                mergeStruct2Fields{idx2} = [intersectFields{idx2},'_s2'];
                        
                warning('Matching struct fields detected; field names of merged struct have been changed.');
            end
                        
        end
        
    end
    
    % Merge myStruct and struct1
    for k = 1:length(struct1Fields);
        
        for j = 1:length(struct1)
            
            myStruct(j).(mergeStruct1Fields{k}) = struct1(j).(struct1Fields{k});
            
        end

    end
    
    for k = 1:length(struct2Fields);
        
        for j = 1:length(struct2)
            
            myStruct(j).(mergeStruct2Fields{k}) = struct2(j).(struct2Fields{k});
            
        end

    end

end
