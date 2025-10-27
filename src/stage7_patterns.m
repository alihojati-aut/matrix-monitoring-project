function result = stage7_patterns(projectRoot, timestamp, params)

    if nargin < 3, params = struct(); end
    if ~isfield(params, 'minSupport'), params.minSupport = 0.01; end
    if ~isfield(params, 'maxLen'), params.maxLen = 3; end
    if ~isfield(params, 'kRec_fd'), params.kRec_fd = 50; end

    reportsDir = fullfile(projectRoot, "reports");
    batchesDir = fullfile(reportsDir, "batches");
    grpDir     = fullfile(reportsDir, "grp");
    ipcaDir    = fullfile(reportsDir, "ipca");
    fdDir      = fullfile(reportsDir, "fd");

    binMatPath = fullfile(reportsDir, sprintf("binary_matrix_%s.mat", timestamp));
    S = load(binMatPath, 'X','txnIds','itemIds');
    X = S.X; itemIds = string(S.itemIds);
    [~, n] = size(X);

    Sgrp = load(fullfile(grpDir, sprintf("grp_model_R_%s.mat", timestamp)), 'R');
    R = Sgrp.R;

    ipcaModelPath = fullfile(ipcaDir, sprintf("ipca_model_final_%s.mat", timestamp));
    if exist(ipcaModelPath, 'file')
        Sip = load(ipcaModelPath, 'ipcaModel');
        if isfield(Sip,'ipcaModel')
            Wt = Sip.ipcaModel.Components;
            mu = Sip.ipcaModel.Mean;
            W = Wt.';
        else
            W = []; mu = [];
        end
    else
        W = []; mu = [];
    end

    X_grp_bin = sparse(0,n);
    X_ipca_bin = sparse(0,n);
    X_fd_bin = sparse(0,n);

    files = dir(fullfile(batchesDir, 'batch_*.mat'));
    for b = 1:numel(files)
        Sb = load(fullfile(batchesDir, files(b).name), 'Xb');
        Xb = Sb.Xb;
        m_b = size(Xb,1);

        Sg = load(fullfile(grpDir, sprintf('grp_batch_%03d.mat', b)), 'Sg');
        Xhat_g = Sg.Sg * R.';
        Xg_bin = binarize_like_original(Xhat_g, Xb);
        X_grp_bin = [X_grp_bin; Xg_bin];

        if ~isempty(W)
            Xb_d = double(full(Xb));
            Xc = bsxfun(@minus, Xb_d, mu);
            Xhat_i = (Xc * (W*W.')) + mu;
        else
            Xb_d = double(full(Xb));
            mu_b = mean(Xb_d,1);
            Xc = bsxfun(@minus, Xb_d, mu_b);
            [~,~,V] = svd(Xc,'econ');
            k = min(50, size(V,2)); Wb = V(:,1:k);
            Xhat_i = (Xc * (Wb*Wb')) + mu_b;
        end
        Xi_bin = binarize_like_original(Xhat_i, Xb);
        X_ipca_bin = [X_ipca_bin; Xi_bin];

        B = load(fullfile(fdDir, sprintf('fd_batch_%03d.mat', b)), 'B'); 
        [~,~,Vb] = svd(double(B.B),'econ');
        kfd = min(params.kRec_fd, size(Vb,2));
        Vk = Vb(:,1:kfd);
        Xhat_f = double(full(Xb)) * (Vk*Vk.');
        Xf_bin = binarize_like_original(Xhat_f, Xb);
        X_fd_bin = [X_fd_bin; Xf_bin];
    end

    minSup = params.minSupport; 
    maxLen = params.maxLen;
    T_orig = fpm_apriori_simple(X,         itemIds, minSup, maxLen);
    T_grp  = fpm_apriori_simple(X_grp_bin, itemIds, minSup, maxLen);
    T_ipca = fpm_apriori_simple(X_ipca_bin,itemIds, minSup, maxLen);
    T_fd   = fpm_apriori_simple(X_fd_bin,  itemIds, minSup, maxLen);

    outDir = fullfile(reportsDir, "patterns");
    if ~exist(outDir,"dir"), mkdir(outDir); end
    p1 = fullfile(outDir, sprintf('patterns_orig_%s.csv', timestamp));
    p2 = fullfile(outDir, sprintf('patterns_grp_%s.csv', timestamp));
    p3 = fullfile(outDir, sprintf('patterns_ipca_%s.csv', timestamp));
    p4 = fullfile(outDir, sprintf('patterns_fd_%s.csv', timestamp));
    writetable(T_orig, p1);
    writetable(T_grp,  p2);
    writetable(T_ipca, p3);
    writetable(T_fd,   p4);

    disp('Stage 7 — Frequent Patterns (preview top 5 of each):');
    disp(head_if(T_orig,5, 'ORIG'));
    disp(head_if(T_grp, 5, 'GRP'));
    disp(head_if(T_ipca,5, 'IPCA'));
    disp(head_if(T_fd,  5, 'FD'));

    result = struct('orig',p1,'grp',p2,'ipca',p3,'fd',p4);
end

function T = head_if(T, k, tag)
    if height(T) > k, T = T(1:k,:); end
    T.method = repmat(string(tag), height(T),1);
end
