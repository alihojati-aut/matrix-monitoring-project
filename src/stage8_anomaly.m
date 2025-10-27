function report = stage8_anomaly(projectRoot, timestamp, b_inject, params)

    if nargin < 4, params = struct(); end
    if ~isfield(params,'nAnom'),       params.nAnom = 20; end
    if ~isfield(params,'patternSize'), params.patternSize = 4; end
    if ~isfield(params,'kRec_fd'),     params.kRec_fd = 50; end

    reportsDir = fullfile(projectRoot, "reports");
    batchesDir = fullfile(reportsDir, "batches");
    grpDir     = fullfile(reportsDir, "grp");
    ipcaDir    = fullfile(reportsDir, "ipca");
    fdDir      = fullfile(reportsDir, "fd");
    outDir     = fullfile(reportsDir, "anomaly");
    if ~exist(outDir,"dir"), mkdir(outDir); end

    Sbin = safe_load(fullfile(reportsDir, sprintf("binary_matrix_%s.mat", timestamp)), {'X','itemIds'});
    Xfull = spones(double(Sbin.X));
    colFreq = full(sum(Xfull, 1));

    Sb = safe_load(fullfile(batchesDir, sprintf("batch_%03d.mat", b_inject)), {'Xb'});
    Xb = Sb.Xb;

    [Xb_mod, injIdx, Xinj] = inject_anomalies_batch(Xb, colFreq, params.nAnom, params.patternSize);
    injPath = fullfile(outDir, sprintf('batch_%03d_injected_%s.mat', b_inject, timestamp));
    save(injPath, 'Xb','Xb_mod','Xinj','injIdx','-v7.3');
    fprintf('\n[Stage 8] Injected %d anomalies into Batch %d (patternSize=%d)\n', ...
        params.nAnom, b_inject, params.patternSize);

    rowResid = @(A,B) sqrt(sum((A - B).^2, 2));

    SgR = safe_load(fullfile(grpDir, sprintf('grp_model_R_%s.mat', timestamp)), {'R'});
    R = SgR.R;
    Xb_mod_d = double(full(Xb_mod));
    Sg_mod = Xb_mod_d * R;
    Xhat_g  = Sg_mod * R.';
    scores_grp = rowResid(Xb_mod_d, Xhat_g);

    Sip = safe_try_load(fullfile(ipcaDir, sprintf('ipca_model_final_%s.mat', timestamp)), {'ipcaModel'});
    if ~isempty(Sip) && isfield(Sip,'ipcaModel') && isstruct(Sip.ipcaModel) && isfield(Sip.ipcaModel,'Components')
        Wt = Sip.ipcaModel.Components;
        mu = Sip.ipcaModel.Mean;
        W  = Wt.';
        Xc = bsxfun(@minus, Xb_mod_d, mu);
        Xhat_i = (Xc * (W*W.')) + mu;
    else

        mu = mean(Xb_mod_d,1);
        Xc = bsxfun(@minus, Xb_mod_d, mu);
        [~,~,V] = svd(Xc,'econ');
        k = min(50, size(V,2));
        W = V(:,1:k);
        Xhat_i = (Xc * (W*W.')) + mu;
    end
    scores_ipca = rowResid(Xb_mod_d, Xhat_i);

    Sbfd = safe_load(fullfile(fdDir, sprintf('fd_batch_%03d.mat', b_inject)), {'B'});
    B = double(Sbfd.B);
    [~,~,Vb] = svd(B,'econ');
    kfd = min(params.kRec_fd, size(Vb,2));
    Vk = Vb(:,1:kfd);
    Xhat_f = Xb_mod_d * (Vk*Vk.');
    scores_fd = rowResid(Xb_mod_d, Xhat_f);

    m_mod = size(Xb_mod_d,1);
    K1 = params.nAnom;
    K2 = min(m_mod, 2*params.nAnom);
    Kp = max(1, round(0.05*m_mod));

    eval_grp = hit_at_k(scores_grp, injIdx, m_mod, [K1,K2,Kp]);
    eval_ipc = hit_at_k(scores_ipca, injIdx, m_mod, [K1,K2,Kp]);
    eval_fd  = hit_at_k(scores_fd,  injIdx, m_mod, [K1,K2,Kp]);

    Tgrp = pack_scores('GRP', scores_grp, injIdx);
    Tipc = pack_scores('IPCA',scores_ipca, injIdx);
    Tfd  = pack_scores('FD',  scores_fd,  injIdx);

    csvG = fullfile(outDir, sprintf('anomaly_scores_grp_b%03d_%s.csv', b_inject, timestamp));
    csvI = fullfile(outDir, sprintf('anomaly_scores_ipca_b%03d_%s.csv', b_inject, timestamp));
    csvF = fullfile(outDir, sprintf('anomaly_scores_fd_b%03d_%s.csv', b_inject, timestamp));
    writetable(Tgrp, csvG); writetable(Tipc, csvI); writetable(Tfd, csvF);

    fprintf('\n[Stage 8] Detection summary on Batch %d (m_mod=%d):\n', b_inject, m_mod);
    pretty_eval('GRP',  eval_grp, K1, K2, Kp);
    pretty_eval('IPCA', eval_ipc, K1, K2, Kp);
    pretty_eval('FD',   eval_fd,  K1, K2, Kp);

    report = struct('scores_grp',csvG,'scores_ipca',csvI,'scores_fd',csvF, ...
                    'batch_injected', injPath);
end


function S = safe_load(matPath, varList, retries, delaySec)
    if nargin<3, retries=12; end
    if nargin<4, delaySec=0.25; end
    for r=1:retries
        try
            S = load(matPath, varList{:});
            return;
        catch ME
            if r==retries
                rethrow(ME);
            else
                pause(delaySec);
            end
        end
    end
end

function S = safe_try_load(matPath, varList, retries, delaySec)
    try
        S = safe_load(matPath, varList, retries, delaySec);
    catch
        S = [];
    end
end

function ev = hit_at_k(scores, injIdx, m, Ks)
    [~, ord] = sort(scores, 'descend');   % ناهنجارتر ⇒ امتیاز بزرگ‌تر
    inj = false(m,1); inj(injIdx) = true;
    ev = struct();
    for i = 1:numel(Ks)
        K = Ks(i);
        topK = ord(1:min(K,m));
        hits = sum(inj(topK));
        ev.(['K' num2str(K)]) = hits;
        ev.(['hitRate_K' num2str(K)]) = hits / max(numel(injIdx),1);
    end
end

function T = pack_scores(method, scores, injIdx)
    m = numel(scores);
    isInj = false(m,1); isInj(injIdx) = true;
    T = table((1:m).', double(scores(:)), logical(isInj), ...
        'VariableNames', {'row','score','injected'});
    T.method = repmat(string(method), m, 1);
    T = T(:, {'method','row','score','injected'});
    T = sortrows(T, 'score', 'descend');
end

function pretty_eval(name, ev, K1, K2, Kp)
    fprintf('  %4s: Hit@%d = %d (%.1f%%) | Hit@%d = %d (%.1f%%) | Hit@%d = %d (%.1f%%)\n', ...
        name, K1, ev.(['K' num2str(K1)]), 100*ev.(['hitRate_K' num2str(K1)]), ...
        K2, ev.(['K' num2str(K2)]), 100*ev.(['hitRate_K' num2str(K2)]), ...
        Kp, ev.(['K' num2str(Kp)]), 100*ev.(['hitRate_K' num2str(Kp)]));
end
