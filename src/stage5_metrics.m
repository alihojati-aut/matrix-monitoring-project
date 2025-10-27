function out = stage5_metrics(cmd, varargin)

persistent allMetrics csvPath matPath

switch lower(string(cmd))
    case "init"
        reportsDir = varargin{1};
        timestamp  = varargin{2};
        csvPath = fullfile(reportsDir, sprintf('metrics_%s.csv', timestamp));
        matPath = fullfile(reportsDir, sprintf('metrics_%s.mat', timestamp));
        allMetrics = table();
        out = struct('csvPath',csvPath,'matPath',matPath);

    case "add"
        method   = varargin{1};
        batchIdx = varargin{2};
        Xb       = varargin{3};
        Xhat     = varargin{4};
        dims     = varargin{5};
        time_s   = varargin{6};

        M = compute_metrics(Xb, 'recon', Xhat);
        row = make_metrics_row(method, batchIdx, Xb, dims, time_s, M);

        if isempty(allMetrics), allMetrics = row; else, allMetrics = [allMetrics; row]; end %#ok<AGROW>
        out = row;

    case "save"
        if isempty(allMetrics)
            warning('No metrics to save. Call stage5_metrics(''add'', ...) first.');
            out = table(); return;
        end
        writetable(allMetrics, csvPath);
        save(matPath, 'allMetrics');
        out = allMetrics;

    otherwise
        error('Unknown command: %s', string(cmd));
end
end
