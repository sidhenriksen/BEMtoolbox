function string = id2string(id);

    id_string = num2str(id);
    K = 1;
    
    N = ceil(length(id_string)/K);
    
    string = [];
    for j = 1:N;
        start = (j-1)*K + 1;
        stop = j*K;
        stop = min([length(id_string),stop]);
        
        current_id = str2double(id_string(start:stop));
        
        current_as_char = char(mod(current_id,122-97)+97);
        string = [string,current_as_char];
    end


end