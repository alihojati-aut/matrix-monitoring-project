function [X, txnIds, itemIds] = build_txn_item_matrix(T)
    assert(all(ismember(["InvoiceNo","StockCode"], string(T.Properties.VariableNames))), ...
        'Columns InvoiceNo and StockCode are required.');

    [txnIds, ~, r]  = unique(T.InvoiceNo, 'stable');   % m×1
    [itemIds, ~, c] = unique(T.StockCode, 'stable');   % n×1

    m = numel(txnIds);
    n = numel(itemIds);

    v = ones(height(T), 1, 'double');
    X = sparse(r, c, v, m, n);

    X = spones(X);
end
