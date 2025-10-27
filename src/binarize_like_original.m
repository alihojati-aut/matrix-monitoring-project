function Xbin = binarize_like_original(Xhat, Xref)
    Xhat = double(full(Xhat));
    if ~issparse(Xref), Xref = sparse(double(Xref)); end

    [m, n] = size(Xhat);

    kPerRow = full(sum(Xref ~= 0, 2));   % m×1
    iiC = cell(m,1);
    jjC = cell(m,1);

    for i = 1:m
        k = kPerRow(i);
        if k <= 0
            iiC{i} = []; jjC{i} = [];
            continue;
        end
        row = Xhat(i, :);

        if exist('maxk','file') == 2
            [~, idx] = maxk(row, min(k, n));
        else
            [~, ord] = sort(row, 'descend');
            idx = ord(1:min(k, n));
        end

        iiC{i} = repmat(i, numel(idx), 1);
        jjC{i} = idx(:);
    end

    if m > 0
        ii = vertcat(iiC{:});
        jj = vertcat(jjC{:});
    else
        ii = []; jj = [];
    end

    v = ones(numel(ii), 1, 'double');
    Xbin = sparse(ii, jj, v, m, n);
    Xbin = spones(Xbin);
end
