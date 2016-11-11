function y = findclosest(x,vec);
    % Usage: y = findclosest(x,vec);
    % Returns the index of vec whose value is closest in to x
    % to x.
    % E.g., findclosest(3,[1.2, 3.2, 4]) would return 2.
    [a,y] = min(abs(vec-x));
end