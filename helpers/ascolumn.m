function y = ascolumn(y);

%    if ~iscolumn(y);
%        y = y';
%    end
    
    if size(y,1) < size(y,2)
        y = y';
    end
    
end