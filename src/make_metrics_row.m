function Trow = make_metrics_row(method, batchIdx, X_batch, dims, time_s, M)

    [m,n] = size(X_batch);
    nnzX  = nnz(X_batch);

    d   = field_or(dims,'d',NaN);
    k   = field_or(dims,'k',NaN);
    ell = field_or(dims,'ell',NaN);
    kRec= field_or(dims,'kRec',NaN);

    Trow = table( ...
        string(method), uint32(batchIdx), uint32(m), uint32(n), uint32(nnzX), ...
        double(d), double(k), double(ell), double(kRec), ...
        double(time_s), ...
        double(M.recErr), double(M.evr), double(M.froX), double(M.froXhat), ...
        'VariableNames', {'method','batch','m','n','nnz','d','k','ell','kRec', ...
                          'time_s','recErr','evr','froX','froXhat'});
end

function v = field_or(S, f, dv)
    if isstruct(S) && isfield(S,f), v = S.(f); else, v = dv; end
end
