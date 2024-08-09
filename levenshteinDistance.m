function distance = levenshteinDistance(str1, str2)
    m = length(str1);
    n = length(str2);
    % Initialize the distance matrix
    D = zeros(m + 1, n + 1);
    for i = 1:m
        D(i + 1, 1) = i;
    end
    for j = 1:n
        D(1, j + 1) = j;
    end
    % Compute the distance
    for j = 1:n
        for i = 1:m
            cost = ~strcmp(str1(i), str2(j));
            D(i + 1, j + 1) = min([D(i, j + 1) + 1, D(i + 1, j) + 1, D(i, j) + cost]);
        end
    end
    % The Levenshtein distance is the value in the bottom-right cell of the matrix
    distance = D(m + 1, n + 1);
end