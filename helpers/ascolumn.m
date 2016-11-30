function y = ascolumn(y);

    if ~iscolumn(y);
        y = y';
    end
    
end