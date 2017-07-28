

function bestK = fit_chi_square(x)

    KMax = 20;

    
    Ps = zeros(1,KMax);
    for k = 1:KMax
            
        z = x/mean(x) * k;
        
        myChi2 = chi2rnd(k,[1,length(x)]);
                
        [~,p] = kstest2(z,myChi2);
        
        Ps(k) = p;
        
        %% This is for computing log likelihood which doesn't really work
        %  when we're doing this funny thing with scaling the means.
        %p = chi2pdf(z,k);
        
        %logLikelihood = sum( log(p) );                
        
        %logLs(k) = logLikelihood;
        
        
    end
    
    [~,bestK] = max(Ps);

end