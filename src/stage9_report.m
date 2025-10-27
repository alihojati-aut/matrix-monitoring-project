function finalT = stage9_report(projectRoot, timestamp, opts)

    if nargin < 3, opts = struct(); end
    if ~isfield(opts,'topN_patterns'), opts.topN_patterns = 50; end
    if ~isfield(opts,'b_inject'),      opts.b_inject      = 0;  end

    reportsDir = fullfile(projectRoot, "reports");
    outDir = fullfile(reportsDir, "final");
    if ~exist(outDir,"dir"), mkdir(outDir); end

    metricsCsv = fullfile(reportsDir, sprintf('metrics_%s.csv', timestamp));
    assert(isfile(metricsCsv), 'Metrics CSV not found: %s', metricsCsv);
    M = readtable(metricsCsv);
    M.method = string(M.method);
    methods = ["GRP","IPCA","FD"];
    recMean = NaN(1,3); evrMean = NaN(1,3); tMean = NaN(1,3);
    for i=1:numel(methods)
        sel = M.method==methods(i);
        recMean(i) = mean(double(M.recErr(sel)),'omitnan');
        evrMean(i) = mean(double(M.evr(sel)),'omitnan');
        tMean(i)   = mean(double(M.time_s(sel)),'omitnan');
    end

    pattDir = fullfile(reportsDir, "patterns");
    pOrig = fullfile(pattDir, sprintf('patterns_orig_%s.csv', timestamp));
    pGRP  = fullfile(pattDir, sprintf('patterns_grp_%s.csv',  timestamp));
    pIPCA = fullfile(pattDir, sprintf('patterns_ipca_%s.csv', timestamp));
    pFD   = fullfile(pattDir, sprintf('patterns_fd_%s.csv',   timestamp));
    pattOverlap = NaN(1,3);
    if isfile(pOrig) && isfile(pGRP) && isfile(pIPCA) && isfile(pFD)
        Torig = readtable(pOrig); Torig.items = string(Torig.items);
        Tgrp  = readtable(pGRP);  Tgrp.items  = string(Tgrp.items);
        Tipca = readtable(pIPCA); Tipca.items = string(Tipca.items);
        Tfd   = readtable(pFD);   Tfd.items   = string(Tfd.items);
        N = min([opts.topN_patterns, height(Torig)]);
        baseSet = Torig.items(1:N);
        pattOverlap(1) = overlap_at_N(Tgrp.items,  baseSet, N);
        pattOverlap(2) = overlap_at_N(Tipca.items, baseSet, N);
        pattOverlap(3) = overlap_at_N(Tfd.items,   baseSet, N);
    end

    hitAtK = NaN(1,3);
    if opts.b_inject > 0
        anomDir = fullfile(reportsDir, "anomaly");
        injMat = fullfile(anomDir, sprintf('batch_%03d_injected_%s.mat', opts.b_inject, timestamp));
        if isfile(injMat)
            Sinj = load(injMat, 'injIdx'); K = numel(Sinj.injIdx);
            fG = fullfile(anomDir, sprintf('anomaly_scores_grp_b%03d_%s.csv',  opts.b_inject, timestamp));
            fI = fullfile(anomDir, sprintf('anomaly_scores_ipca_b%03d_%s.csv', opts.b_inject, timestamp));
            fF = fullfile(anomDir, sprintf('anomaly_scores_fd_b%03d_%s.csv',   opts.b_inject, timestamp));
            if isfile(fG), hitAtK(1) = compute_hit_atK(fG, Sinj.injIdx, K); end
            if isfile(fI), hitAtK(2) = compute_hit_atK(fI, Sinj.injIdx, K); end
            if isfile(fF), hitAtK(3) = compute_hit_atK(fF, Sinj.injIdx, K); end
        end
    end

    finalT = table(methods.', recMean.', evrMean.', tMean.', pattOverlap.', hitAtK.', ...
        'VariableNames', {'method','recErr_mean','evr_mean','time_mean','pattern_overlap_atN','hit_atK'});

    R = zeros(height(finalT),1);
    R = R + rank_pos(finalT.recErr_mean, 'ascend');
    R = R + rank_pos(finalT.time_mean,   'ascend');
    R = R + rank_pos(finalT.evr_mean,    'descend');
    if all(~isnan(finalT.pattern_overlap_atN)), R = R + rank_pos(finalT.pattern_overlap_atN, 'descend'); end
    if all(~isnan(finalT.hit_atK)),            R = R + rank_pos(finalT.hit_atK, 'descend');             end
    finalT.overall_rank = R;

    outCsv = fullfile(outDir, sprintf('final_report_%s.csv', timestamp));
    outMat = fullfile(outDir, sprintf('final_report_%s.mat', timestamp));
    writetable(finalT, outCsv);
    save(outMat, 'finalT');
    fprintf('\n[Stage 9] Final report saved:\n%s\n', outCsv);
    disp('Summary (by method):'); disp(finalT);
end

function ov = overlap_at_N(items, baseSet, N)
    items = items(1:min(N,numel(items)));
    ov = sum(ismember(items, baseSet)) / max(N,1);
end
function hit = compute_hit_atK(scoreCsv, injIdx, K)
    T = readtable(scoreCsv);
    T = sortrows(T, 'score', 'descend');
    topK = T.row(1:min(K,height(T)));
    inj = false(max([max(T.row), max(injIdx)]),1);
    inj(injIdx) = true;
    hit = sum(inj(topK)) / max(numel(injIdx),1);
end
function r = rank_pos(x, direction)
    x = x(:); isN = isnan(x);
    switch direction
        case 'ascend',  [~,ord] = sort(x,'ascend','MissingPlacement','last');
        case 'descend', [~,ord] = sort(x,'descend','MissingPlacement','last');
    end
    r = zeros(size(x)); r(ord) = 1:numel(x);
    r(isN) = numel(x);
end
