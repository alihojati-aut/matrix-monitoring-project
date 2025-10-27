function T = load_preprocess(csvPath)

    assert(isfile(csvPath), 'File not found: %s', csvPath);
    opts = detectImportOptions(csvPath, 'TextType','string');
    T = readtable(csvPath, opts);

    mustHave = ["InvoiceNo","StockCode","InvoiceDate"];
    vn = string(T.Properties.VariableNames);
    missingCols = setdiff(mustHave, intersect(mustHave, vn));
    assert(isempty(missingCols), 'Missing required columns: %s', strjoin(missingCols, ', '));

    if ~isstring(T.InvoiceNo), T.InvoiceNo = string(T.InvoiceNo); end
    if ~isstring(T.StockCode), T.StockCode = string(T.StockCode); end
    if ~isdatetime(T.InvoiceDate)

        try
            T.InvoiceDate = datetime(T.InvoiceDate, 'InputFormat','dd/MM/yyyy HH:mm');
        catch
            T.InvoiceDate = datetime(T.InvoiceDate);
        end
    end

    keyOk = ~ismissing(T.InvoiceNo) & ~ismissing(T.StockCode) & ~ismissing(T.InvoiceDate);
    T = T(keyOk, :);

    isReturn = startsWith(T.InvoiceNo, "C", 'IgnoreCase', true);
    T = T(~isReturn, :);

    if any(strcmpi('Quantity', T.Properties.VariableNames))
        T = T(T.Quantity > 0, :);
    end

end
