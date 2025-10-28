function batches = stream_batches(X, batchSize)
    m = size(X,1);
    if batchSize <= 0
        error('batchSize must be a positive integer.');
    end

    starts = 1:batchSize:m;
    nB = numel(starts);
    batches = cell(nB,1);

    for i = 1:nB
        r1 = starts(i);
        r2 = min(r1 + batchSize - 1, m);
        batches{i} = X(r1:r2, :);
    end
end
