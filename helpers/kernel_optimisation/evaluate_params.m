function [cost,binocularKernel] = evaluate_params(NimModel,params,is2d,costType)
    % EVALUATE_PARAMS   Wrapper for contrast_cost.
    % Evaluates the current Gabor parameters and returns
    % mean spike count difference for high/low contrast.
    %
    % Parameters
    % ----------
    % NimModel : NimModel struct
    % params : (Binocular) Gabor parameters to evaluate
    %
    % Returns
    % --------
    % cost : the mean spike count difference between high and low contrast
    % frames, as determined by the parameterised Gabor
    % binocularGabor : the binocular Gabor
    
    if nargin < 3
        is2d = false;
    end
    
    if nargin < 4
        costType = 2;
    end

    x = 1:(size(NimModel.stim,2)/2);
    
    if is2d 
        t = 1:NimModel.NimFit.stim_params.dims(1);
        
        binocularKernel = binocular_spacetime_gaussian(x,t,params);
    else
        binocularKernel = binocular_gabor(x,params);
    end
    
    
    cost = contrast_cost(NimModel,binocularKernel,costType);
    
end

function monocularGabor = monocular_gabor(x,params)

    A = params(1); % amplitude
    m = params(2); % mean
    s = params(3); % sd
    f = params(4); % frequency
    phi = params(5); % phase
    
    monocularGabor = A .* exp(-(x-m).^2 / (2*s.^2)) .* cos(2*pi*f*(x-m) + phi);


end

function binocularGabor = binocular_gabor(x,params)

    leftParams = params(1:5);
    
    rightParams = params(6:end);
    
    leftGabor = monocular_gabor(x,leftParams);
    
    rightGabor = monocular_gabor(x,rightParams);
    
    binocularGabor = [leftGabor,rightGabor]; 
end


function spacetimeGaussian = spacetime_gaussian(x,t,params)

    A = params(1); % amplitude
    mx = params(2); % mean (space)
    sx = params(3); % sd (space)
    mt = params(4); % mean (time)
    st = params(5); % sd (time)
    r = params(6); % correlation
    
    
    [X,T] = meshgrid(x,t);
        
    spacetimeGaussian = A * exp(- 1/(2*(1-r^2)) *( ...
        (X(:)-mx).^2 ./ sx^2 + ...
        (T(:)-mt).^2 ./ st^2 + ...
        2*r*(X(:)-mx).*(T(:)-mt)./(sx*st)));
        
    
end

function binocularSpacetimeGaussian = binocular_spacetime_gaussian(x,t,params)


    leftParams = params(1:6);
    
    rightParams = params(7:end);
    
    leftSpacetimeGaussian = spacetime_gaussian(x,t,leftParams);
    
    rightSpacetimeGaussian = spacetime_gaussian(x,t,rightParams);
    
    binocularSpacetimeGaussian = [leftSpacetimeGaussian;rightSpacetimeGaussian]';

end