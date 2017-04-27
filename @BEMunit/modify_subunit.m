function varargout = modify_subunit(bem,k_sub,varargin)
    % Usage: bem = bem.modify_subunit(k_sub,<attributes>)
    % Where k_sub is the index of the subunit to be modified
    %
    % ===Example usage===
    % Change left and right RF phase:
    % bem = bem.modify_subunit(2,'L_phi',pi/2,'R_phi',0);
    % 
    % Change output NL to half-squaring (here Pos would be a
    % pre-defined half-wave rectifying function):
    % bem = bem.modify_subunit(1,'NL',@(x)(Pos(x).^2));

    rf = bem.subunits(k_sub).rf_params;

    for j = 1:2:(length(varargin)-1);
        if strcmp('L_',varargin{j}(1:2))
            assert(isfield(rf.left,varargin{j}(3:end)),sprintf('Invalid RF parameter given: %s\n',varargin{j}));
            rf.left.(varargin{j}(3:end)) = varargin{j+1};
        elseif strcmp('R_',varargin{j}(1:2));
            assert(isfield(rf.right,varargin{j}(3:end)),sprintf('Invalid RF parameter given: %s\n',varargin{j}));
            rf.right.(varargin{j}(3:end)) = varargin{j+1};
        elseif strcmpi('NL',varargin{j});
            rf.NL = varargin{j+1};                    
        else
            if ~isfield(rf,varargin{j});
                assert(isfield(rf.left,varargin{j}),sprintf('Invalid RF parameter given: %s\n',varargin{j}));
            end
            rf.(varargin{j}) = varargin{j+1};
        end
    end

    bem.subunits(k_sub).rf_params = rf;
    bem.update();
    
    if nargout
        varargout = {bem};
    end
    
end