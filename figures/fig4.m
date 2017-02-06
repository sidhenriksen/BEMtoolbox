

function fig4()
    constants = figure_constants();

    allCells = get_plotter_data();
    
    slopes = arrayfun(@(K)get_slope(K),allCells);
    rs = arrayfun(@(K)get_r(K),allCells);
    
    figure(); 
    plot(rs,slopes,'k o','markersize',constants.markersize)

end

function slope = get_slope(myCell)


    externalVariance = strip_uc(myCell.totalVariance-myCell.internalVariance);
    
    spikeCount = strip_uc(myCell.twopassSpikeCount);
    
    externalFanoFactor = externalVariance./spikeCount;
    
    [~,slope] = regression2(externalFanoFactor,spikeCount);
    
end

function r = get_r(myCell)

    externalVariance = strip_uc(myCell.totalVariance-myCell.internalVariance);
    
    spikeCount = strip_uc(myCell.twopassSpikeCount);
    
    externalFanoFactor = externalVariance./spikeCount;
    
    r = regression2(externalFanoFactor,spikeCount);

end