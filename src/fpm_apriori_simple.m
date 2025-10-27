function outT = fpm_apriori_simple(Xbin, itemIds, minSupport, maxLen)
    if nargin < 3 || isempty(minSupport), minSupport = 0.01; end
    if nargin < 4 || isempty(maxLen),     maxLen     = 3;    end
    maxLen = min(max(maxLen,1),3);

    if iscell(itemIds), itemIds = string(itemIds); end
    itemIds = string(itemIds);
    Xbin = spones(Xbin) ~= 0;

    [m, n] = size(Xbin);
    minCount = ceil(minSupport * m);

    outT = table();

    c1 = full(sum(Xbin,1));
    keep1 = find(c1 >= minCount);
    if isempty(keep1)
        if maxLen == 1
            outT = empty_table();
        end
        return;
    end
    L1 = keep1(:);
    C1 = c1(L1).';

    T1 = pack_rows(1, L1, C1, itemIds, m);
    outT = [outT; T1]; %#ok<AGROW>
    if maxLen == 1, return; end

    if numel(L1) >= 2
        pairs = nchoosek(L1, 2);
        cnt2  = count_pairs(Xbin, pairs);
        keep2 = find(cnt2 >= minCount);
        L2 = pairs(keep2, :);
        C2 = cnt2(keep2);
        if ~isempty(L2)
            T2 = pack_rows(2, L2, C2, itemIds, m);
            outT = [outT; T2]; %#ok<AGROW>
        end
    else
        if maxLen >= 2
            return;
        end
    end
    if maxLen == 2, return; end
    if exist('L2','var') && ~isempty(L2)
        triples = gen_triples_from_L2(L2);
        if ~isempty(triples)
            cnt3  = count_triples(Xbin, triples);
            keep3 = find(cnt3 >= minCount);
            L3 = triples(keep3, :);
            C3 = cnt3(keep3);
            if ~isempty(L3)
                T3 = pack_rows(3, L3, C3, itemIds, m);
                outT = [outT; T3]; %#ok<AGROW>
            end
        end
    end

    if ~isempty(outT)
        outT = sortrows(outT, {'len','support','count'}, {'ascend','descend','descend'});
    else
        outT = empty_table();
    end
end


function cnt = count_pairs(X, pairs)
    k = size(pairs,1);
    cnt = zeros(k,1);
    for i=1:k
        a = pairs(i,1); b = pairs(i,2);
        cnt(i) = full(sum(X(:,a) & X(:,b)));
    end
end

function triples = gen_triples_from_L2(L2)
    if isempty(L2), triples = zeros(0,3); return; end
    L2 = sort(L2,2);
    triples = [];
    for i=1:size(L2,1)
        for j=i+1:size(L2,1)
            a1=L2(i,1); b1=L2(i,2);
            a2=L2(j,1); b2=L2(j,2);
            if a1==a2 && b1~=b2
                t = sort([a1, b1, b2]);
                triples = [triples; t]; %#ok<AGROW>
            end
        end
    end
    if isempty(triples)
        triples = zeros(0,3);
    else
        triples = unique(triples, 'rows');
    end
end

function cnt = count_triples(X, triples)
    k = size(triples,1);
    cnt = zeros(k,1);
    for i=1:k
        a=triples(i,1); b=triples(i,2); c=triples(i,3);
        cnt(i) = full(sum(X(:,a) & X(:,b) & X(:,c)));
    end
end

function T = pack_rows(len, idxMat, counts, itemIds, m)
    if isempty(idxMat)
        T = empty_table(); return;
    end
    k = size(idxMat,1);

    counts = double(counts(:));
    if numel(counts) ~= k
        error('pack_rows: counts length (%d) does not match number of patterns (%d).', numel(counts), k);
    end

    itemsStr = strings(k,1);
    switch len
        case 1
            itemsStr = itemIds(idxMat);
        otherwise
            for i=1:k
                itemsStr(i) = strjoin(itemIds(idxMat(i,:)), ", ");
            end
    end

    support = counts / max(m,1);
    T = table( ...
        repmat(uint8(len), k, 1), ...
        itemsStr, ...
        support, ...
        counts, ...
        'VariableNames', {'len','items','support','count'});
end

function T = empty_table()
    T = table( uint8([]), string([]), double([]), double([]), ...
        'VariableNames', {'len','items','support','count'});
end
