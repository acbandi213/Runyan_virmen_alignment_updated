function bestTrialIndices = findBestMatchingTrials(trueConditions,estimatedConditions,startTrial)
% Convert cell arrays to strings for comparison
    trueStr = strjoin(trueConditions, '|'); % '|' as a delimiter
    estimatedStr = strjoin(estimatedConditions, '|'); % '|' as a delimiter
    % Get individual trials by splitting the strings
    trueTrials = strsplit(trueStr, '|');
    estimatedTrials = strsplit(estimatedStr, '|');
    % Define the window size
    windowSize = length(estimatedTrials);
    % Initialize variables
    numTrueTrials = length(trueTrials);
    bestDistance = Inf;
    bestTrialIndices = [];
    % Iterate through true trials with a sliding window
    for i = startTrial:numTrueTrials - windowSize + 1
        trueWindow = trueTrials(i:i+windowSize-1);
        % Calculate Levenshtein distance for the window against the estimated trials
        windowDistance = 0;
        for j = 1:windowSize
            windowDistance = windowDistance + levenshteinDistance(trueWindow{j}, estimatedTrials{j});
        end
        % Update best match if the distance is smaller
        if windowDistance < bestDistance
            bestDistance = windowDistance;
            bestTrialIndices = [i, i + windowSize - 1]; % Store indices of best matching window
        end
    end
    fprintf('Best matching true trials window indices: %d to %d\n', bestTrialIndices(1), bestTrialIndices(2));
    fprintf('Total Levenshtein distance for best-matching window: %d\n', bestDistance);
    end