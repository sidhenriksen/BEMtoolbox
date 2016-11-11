function y = strip_uc(x);
    y = [x([1,3],:)]; y = y(:);
    y = [y;x(2,1)];
end
