function [L,R] = make_rls(K,N,dx,bc)
    % make_rls  creates rls
    % Parameters
    % K : integer (number of rows)
    % N : integer (number of columns)
    % dx : disparity (in pixels)
    % bc: binocular correlation(-1, 0, or 1)

    assert(round(dx)==dx,'Error: disparity must be integer.')
    
    assert(ismember(bc,[-1,0,1]),'Error: invalid correlation value.')
            
    L = randi(3,[K,N])-2;
    
    if bc == 0
        
        R = randi(3,[K,N])-2;
        
    else
        
        R = latshift3(L*bc,dx);
        
    end
    
    
end