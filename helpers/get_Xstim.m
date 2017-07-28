

function Xstim = get_Xstim(NimModel)

    stim_params = NimModel.NimFit.stim_params;
    
    Xstim = NIM.create_time_embedding(NimModel.stim,stim_params);

end