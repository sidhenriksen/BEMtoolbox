function [NimModel,LLCv,LLTrain] = fit_NimModel(NimCell,nExc,nInh,optStruct)
    % fit_NImModel Fit a NIM model using a specifiec architecture.
    % Usage:
    % [NimModel,LL] = fit_NimModel(NimCell,n_exc,n_inh,optStruct)
    % 
    % NimCell : a NimCell structure (obtained from make_NimCell)
    % n_exc : number of excitatory subunits
    % n_inh : number of inhibitory/suppressive subunits
    % runDur : run duration (e.g. 10, 30, 750 ms)
    % optStruct : 
    %
    % Returns a NimModel and the cross-validated log likelihood   
     
    

    % accepted is rectpow, quad, rectlin, lin, nonpar
    nlType = 'rectpow';
    
    spkNL = 'softplus';
    
    runDur = 30;
    
    nLags = 15;
    
    propData = 1.0; % proportion of data to use
    
    pTrain = 1.0;
        
    silent = 0;
        
    smoothLambda = 500; %spacetime smoothness regularisation
        
    L1Lambda=1; % L1 norm regularisation

    
    if nargin == 4
        % this can be used to override the previous parameters
        unpackStruct(optStruct);
    end    
    

    

    dt = NimCell.times(2)-NimCell.times(1);
    
    stim = NimCell.stim; 
    spiketimes = NimCell.spiketimes;
    times = NimCell.times; 
    duration = NimCell.duration;
    
    % mask for which "frames" to include.. we exclude the first 200 ms
    % of every trial for a variety of reasons    
    include = logical((times >= 0.2) .* (duration == runDur) .* (abs(NimCell.dxs)<1e2));
    
    % if propData < 1, then only use a fraction of the data
    if propData < 1
        prop_mask = rand([1,length(include)]) < propData;
        include = logical(include.*prop_mask);
    end
    
    if ~isfield(optStruct,'indexTrain')        
        [indexTrain,indexCv] = partition_data(times(include),pTrain);
    else
        logicalIndex = zeros(1,length(include));
        logicalIndex(optStruct.indexTrain) = 1;
        logicalIndex(optStruct.indexCv) = 2;
        logicalIndex = logicalIndex(include);
        indexTrain = find(logicalIndex==1);
        indexCv = find(logicalIndex==2);
        
    end
    
    nPix = size(stim,2);
    
    params_stim = NIM.create_stim_params([nLags,nPix],'stim_dt',dt);

    XcellAll{1} = NIM.create_time_embedding(stim,params_stim);    
    Xcell{1} = XcellAll{1}(include,:); % get rid of the data we don't want to include
    
    rObsAll = histc(spiketimes,(0:(length(stim)-1))*params_stim.dt);    
    rObs = rObsAll(include);
    

    mod_signs = [ones(1,nExc),ones(1,nInh)*-1];

    
    if strcmp(nlType,'quad')
        % we need to add a linear subunit if we're using a quadratic model
        %nlTypes = ['lin',repmat({'quad'},[1,nExc-1 + nInh])];
        
        nlTypes = repmat({nlType},[1,nExc + nInh]);
        
        NimFit = NIM(params_stim,nlTypes,mod_signs,'spkNL',spkNL);            
        NimFit = NimFit.set_reg_params('d2xt',smoothLambda,'l1',L1Lambda);
        NimFit.stim_params.split_pts = [2, nPix/2, 0]; % sets boundary on 2nd dimension, at pixel nPix/2 (middle) and smooths to 0 at this boundary    

        NimFit = NimFit.fit_filters(rObs,Xcell,indexTrain);
        NimFit.fit_history = [];
        NimFit = NimFit.fit_spkNL(rObs,Xcell,indexTrain);


    elseif any(strcmp(nlType,{'rectpow','rectlin','lin','nonpar'}))
        % same as previous, just without adding the linear subunit
        
        nlTypes = repmat({nlType},[1,nExc + nInh]); % sets subunit NL types
        
        offsets = ones(size(nlTypes)); % toggles threshold fitting
 
        NimFit = NIM(params_stim,nlTypes,mod_signs,'spkNL',spkNL); % 
        
        NimFit = NimFit.set_reg_params('d2xt',smoothLambda,'l1',L1Lambda);
                        
        NimFit.stim_params.split_pts = [2, nPix/2, 0]; % sets boundary on 2nd dimension, at pixel nPix/2 (middle) and smooths to 0 at this boundary    


        NimFit = NimFit.fit_filters(rObs,Xcell,indexTrain,'silent',silent,'fit_offsets',offsets);        
        NimFit.fit_history=[];

        NimFit = NimFit.fit_spkNL(rObs,Xcell,indexTrain,'silent',silent);
    else
        error('Unrecognised nonlinearity.');
    end

    [~,rEst] = NimFit.eval_model(rObsAll,XcellAll);    
    LLTrain = NimFit.eval_model(rObs,Xcell,indexTrain);
    LLCv = NimFit.eval_model(rObs,Xcell,indexCv);
    
    NimModel = NimCell;
    NimModel.NimFit=NimFit;
    NimModel.spiketimes = [];
    NimModel.rEst = rEst;
    NimModel.runDur = runDur;
    NimModel.indexTrain = indexTrain;
    NimModel.indexCv = indexCv;

end