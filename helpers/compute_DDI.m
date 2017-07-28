function DDI = compute_DDI(unitData)
            
    totalVariance = strip_uc(unitData.totalSqrtVariance);
    
    spikeCount = strip_uc(unitData.twopassSqrtSpikeCount);
    
    rMax = max(spikeCount);
    
    rMin = min(spikeCount);
    
    RMSerror = sqrt(mean(totalVariance));
    
    DDI = (rMax-rMin)/(rMax-rMin + 2*RMSerror);
    
    
end
