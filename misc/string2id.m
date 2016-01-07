function id = string2id(string);

    
    K = 1;
    
    N = ceil(length(string)/K);
    
    id = [];
    for j = 1:N;
        start = (j-1)*K + 1;
        stop = j*K;
        stop = min([length(string),stop]);
        
        current_id = double(string(start:stop));                
        
        current_as_char = num2str(current_id-97);
        id = [id,current_as_char];
    end


end