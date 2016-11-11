function [indexTrain,indexCv] = partition_data(times,pTrain)
    % Partition data based on trials not frames. This is important since
    % the frames are not actually independent.
    % 
    % Usage: [indexTrain,indexCv] = partition_data(times,pTrain);
    %
    % times : array containing the trial time for each frame relative to
    % the start of the trial
    % 
    % pTrain : proportion of data that goes into the training set

    if nargin < 2
        pTrain = 0.75;
    end

    
    % Assign each frame to its corresponding trial
    framesTrials = cumsum(diff(times) < 0 )+1;
    framesTrials = [1,framesTrials]; % Prepend a 1 to account for diff making the vector of length n-1

    % Split into training and cv set; we do this on a trial-by-trial basis
    nTrials = length(unique(framesTrials));
    
    nTrain = round(nTrials*pTrain);
    rand_index = randperm(nTrials);
    trial_index_train = rand_index(1:nTrain);
    trial_index_cv = rand_index((nTrain+1):end);

    indexTrain = [];
    indexCv = [];
    k = 0;
    for j = 1:nTrials
        add_index = find(framesTrials == j);
        k = k+length(add_index);
        if any(j == trial_index_train);
            indexTrain = [indexTrain,add_index];

        elseif any(j == trial_index_cv);
            indexCv = [indexCv,add_index];

        end
    end
    indexTrain = sort(indexTrain);
    indexCv = sort(indexCv);
end