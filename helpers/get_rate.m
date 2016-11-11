function rObs = get_rate(NimCell)

    dt = NimCell.times(2)-NimCell.times(1);
    
    rObs = histc(NimCell.spiketimes,(0:(length(NimCell.stim)-1))*dt);
end