function [Xb_mod, injIdx, Xinj] = inject_anomalies_batch(Xb, colFreqGlobal, nAnom, patternSize)
    Xb = spones(double(Xb));
    [mb, n] = size(Xb);

    if iscolumn(colFreqGlobal), colFreqGlobal = colFreqGlobal.'; end
    if numel(colFreqGlobal) ~= n
        error('colFreqGlobal must have length n (number of columns in X).');
    end
    colFreqGlobal = double(colFreqGlobal);

    nRare = max(10, round(0.10 * n));
    [~, order] = sort(colFreqGlobal, 'ascend');
    rareCols = order(1:nRare);

    ii = []; jj = [];
    for r = 1:nAnom
        pick = randsample(rareCols, min(patternSize, numel(rareCols)), false);
        ii = [ii; repmat(r, numel(pick), 1)]; %#ok<AGROW>
        jj = [jj; pick(:)]; %#ok<AGROW>
    end
    Xinj = sparse(ii, jj, 1.0, nAnom, n);   % nAnom × n

    Xb_mod = [Xb; Xinj];
    injIdx = (mb+1 : mb+nAnom).';
end
