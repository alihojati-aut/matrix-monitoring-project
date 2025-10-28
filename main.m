function main()
clc; close all; format compact;

%% ====================== تنظیم مسیر پروژه و پارامترها ======================
projectRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(projectRoot, "src"));

dataPath   = fullfile(projectRoot, "data", "OnlineRetail.xlsx");
reportsDir = fullfile(projectRoot, "reports");
dirs = {reportsDir, fullfile(reportsDir,"batches"), fullfile(reportsDir,"grp"), ...
        fullfile(reportsDir,"ipca"), fullfile(reportsDir,"fd"), ...
        fullfile(reportsDir,"figures"), fullfile(reportsDir,"patterns"), ...
        fullfile(reportsDir,"anomaly"), fullfile(reportsDir,"final")};
for i=1:numel(dirs), if ~exist(dirs{i},"dir"), mkdir(dirs{i}); end, end

timestamp = datestr(now,'yyyymmdd_HHMMSS');

batchSize = 1000;
d_grp   = 300;  rngSeed = 42;
k_ipca  = 250;
ell_fd  = 120;  kRec_fd = 50;
minSupport = 0.02;  maxLen = 3;
b_inject = 7;
nAnom = 20;  patternSize = 4;

%% =============== مرحله 1: بارگذاری و پاکسازی داده (Table تمیز) ===============
disp("================================================================");
disp("مرحله 1: بارگذاری و پاک‌سازی داده‌ها");
T = load_preprocess(dataPath);

cleanMatPath = fullfile(reportsDir, sprintf('cleaned_table_%s.mat', timestamp));
save(cleanMatPath, 'T', '-v7.3');
fprintf('Saved cleaned table → %s\n\n', cleanMatPath);

disp('— پیش‌نمایش 5 ردیف از جدول تمیزشده —');
try
    disp(T(1:min(5,height(T)), :));
catch
    head(T);
end

%% ======== مرحله 2: ساخت ماتریس باینری تراکنش–آیتم (X، txnIds، itemIds) ========
disp("================================================================");
disp("مرحله 2: ساخت ماتریس تراکنش–آیتم (باینری)");
[X, txnIds, itemIds] = build_txn_item_matrix(T);
[m, n] = size(X); nz = nnz(X);
fprintf('X: %d×%d  (nonzeros=%d, density=%.6f)\n', m,n,nz,nz/(m*n));
binMatPath = fullfile(reportsDir, sprintf('binary_matrix_%s.mat', timestamp));
save(binMatPath, 'X', 'txnIds', 'itemIds', '-v7.3');
fprintf('Saved binary matrix → %s\n\n', binMatPath);

disp('— پیش‌نمایش 5×10 از —');
disp(full(X(1:min(5,m), 1:min(10,n))));

%% ================= مرحله ۳: تقسیم به بچ‌ها (Streaming Batches) =================
disp("================================================================");
disp("مرحله 3: تقسیم به بچ‌ها");
batches = stream_batches(X, batchSize);
nB = numel(batches);
fprintf('تعداد بچ‌ها = %d  (batchSize = %d)\n', nB, batchSize);
for b = 1:nB
    Xb = batches{b};
    [mb, nb] = size(Xb); nbz = nnz(Xb);
    fprintf('  • Batch %d → %d×%d (nonzeros=%d, density=%.6f)\n', ...
        b, mb, nb, nbz, nbz/(mb*nb));
    save(fullfile(reportsDir,"batches",sprintf('batch_%03d.mat', b)), 'Xb','-v7.3');
end
if nB>0
    disp('— پیش‌نمایش 5×10 از Batch #1 —');
    Xb1 = batches{1};
    disp(full(Xb1(1:min(5,size(Xb1,1)), 1:min(10,size(Xb1,2)))));
end

%% =================== مرحله ۴: اجرای سه روش GRP / IPCA / FD ====================
disp("================================================================");
disp("مرحله 4: اجرای GRP, IPCA, FD روی بچ‌ها و ذخیره خروجی‌ها");

grpTime  = nan(nB,1);
ipcaTime = nan(nB,1);
fdTime   = nan(nB,1);

% --- GRP ---
grpDir = fullfile(reportsDir,"grp");
fprintf('\n=== [GRP] d=%d ===\n', d_grp);
for b = 1:nB
    Xb = batches{b};
    t0 = tic;
    [Sg, Xhat_g, R] = sketch_grp(Xb, d_grp, rngSeed);
    dt = toc(t0);
    grpTime(b) = dt;
    recErr = norm(double(full(Xb)) - Xhat_g, 'fro') / max(norm(double(full(Xb)),'fro'), eps);
    fprintf('  • GRP Batch %d: Sg=%d×%d, err=%.4f, time=%.3fs\n', ...
        b, size(Sg,1), size(Sg,2), recErr, dt);
    save(fullfile(grpDir, sprintf('grp_batch_%03d.mat', b)), 'Sg','R','Xhat_g','-v7.3');

    if b==1
        disp('پیش‌نمایش 5×10 از Sg (Batch #1):');
        disp(Sg(1:min(5,size(Sg,1)), 1:min(10,size(Sg,2))));
    end
end
save(fullfile(grpDir, sprintf('grp_model_R_%s.mat', timestamp)), 'R','-v7.3');
fprintf('[GRP] خروجی‌ها در: %s\n', grpDir);

% --- IPCA ---
ipcaDir = fullfile(reportsDir,"ipca");
fprintf('\n=== [IPCA] k=%d ===\n', k_ipca);
ipcaModel = [];
for b = 1:nB
    Xb = batches{b};
    t0 = tic;
    [Zi, Xhat_i, ipcaModel] = sketch_ipca(Xb, k_ipca, ipcaModel);
    dt = toc(t0);
    ipcaTime(b) = dt;
    recErr = norm(double(full(Xb)) - Xhat_i, 'fro') / max(norm(double(full(Xb)),'fro'), eps);
    fprintf('  • IPCA Batch %d: Zi=%d×%d, err=%.4f, time=%.3fs\n', ...
        b, size(Zi,1), size(Zi,2), recErr, dt);
    save(fullfile(ipcaDir, sprintf('ipca_batch_%03d.mat', b)), 'Zi','-v7.3');
    if b==1
        disp('پیش‌نمایش 5×10 از Zi (Batch #1):');
        disp(Zi(1:min(5,size(Zi,1)), 1:min(10,size(Zi,2))));
    end
end
save(fullfile(ipcaDir, sprintf('ipca_model_final_%s.mat', timestamp)), 'ipcaModel','-v7.3');
fprintf('[IPCA] خروجی‌ها در: %s\n', ipcaDir);

% --- FD ---
fdDir = fullfile(reportsDir,"fd");
fprintf('\n=== [FD] ell=%d (kRec=%d) ===\n', ell_fd, kRec_fd);
for b = 1:nB
    Xb = batches{b};
    t0 = tic;
    [B, Xhat_f] = sketch_fd(Xb, ell_fd, kRec_fd);
    dt = toc(t0);
    fdTime(b) = dt;
    recErr = norm(double(full(Xb)) - Xhat_f, 'fro') / max(norm(double(full(Xb)),'fro'), eps);
    fprintf('  • FD  Batch %d: B=%d×%d, err=%.4f, time=%.3fs\n', ...
        b, size(B,1), size(B,2), recErr, dt);
    save(fullfile(fdDir, sprintf('fd_batch_%03d.mat', b)), 'B','-v7.3');
    if b==1
        disp('پیش‌نمایش 5×10 از B (Batch #1):');
        disp(B(1:min(5,size(B,1)), 1:min(10,size(B,2))));
    end
end
fprintf('[FD] خروجی‌ها در: %s\n', fdDir);

%% ==================== مرحله ۵: محاسبه و ذخیرهٔ شاخص‌ها (Metrics) ====================
disp("================================================================");
disp("مرحله 5: محاسبه و ذخیره شاخص‌ها (همه روش‌ها)");

stage5_metrics('init', reportsDir, timestamp);
clear R
% GRP metrics

for b = 1:nB
    Xb = batches{b};
    S = load(fullfile(grpDir, sprintf('grp_batch_%03d.mat', b)), 'Sg','R','Xhat_g');

    if isfield(S, 'Xhat_g')
        Xhat_g = S.Xhat_g;
    else
        Xhat_g = S.Sg * S.R.';
    end

    stage5_metrics('add','GRP', b, Xb, Xhat_g, struct('d', d_grp), grpTime(b));
end


% IPCA metrics
for b = 1:nB
    Xb   = batches{b};
    Xb_d = double(full(Xb));
    if exist('ipcaModel','var') && isstruct(ipcaModel) && isfield(ipcaModel,'Components')
        Wt = ipcaModel.Components;         % (k × n)
        mu = ipcaModel.Mean;               % (1 × n)
        W  = Wt.';                         % (n × k)
        Xhat_i = (bsxfun(@minus, Xb_d, mu)) * (W*W.') + mu;
    else
        mu_b = mean(Xb_d,1);
        Xc   = bsxfun(@minus, Xb_d, mu_b);
        [~,~,V] = svd(Xc,'econ');
        k = min(k_ipca, size(V,2));
        Wb = V(:,1:k);
        Xhat_i = (Xc * (Wb*Wb')) + mu_b;
    end
    stage5_metrics('add','IPCA',b,Xb,Xhat_i,struct('k',k_ipca), ipcaTime(b));
end

% FD metrics
for b = 1:nB
    Xb = batches{b};
    B = load(fullfile(fdDir, sprintf('fd_batch_%03d.mat', b)), 'B'); 
    B = double(B.B);
    [~,~,Vb] = svd(B,'econ'); 
    kfd = min(kRec_fd, size(Vb,2)); 
    Vk = Vb(:,1:kfd);
    Xhat_f = double(full(Xb)) * (Vk*Vk.');
    stage5_metrics('add','FD',b,Xb,Xhat_f,struct('ell',ell_fd,'kRec',kRec_fd), fdTime(b));
end

All = stage5_metrics('save');
disp('— پیش‌نمایش 5 ردیف از جدول متریک‌ها —');
disp(All(1:min(5,height(All)), :));

%% ==================== مرحله ۶: ترسیم نمودارها (درون main) =======================
disp("================================================================");
disp("مرحله 6: ترسیم نمودارها (درون main)");

T = All;
T.method = string(T.method);
Tbatch = double(T.batch);
methods = unique(T.method,'stable');
figDir = fullfile(reportsDir,'figures'); if ~exist(figDir,'dir'), mkdir(figDir); end
titleSuffix = datestr(now,'yyyy-mm-dd');

if all(ismember({'froX','froXhat'}, T.Properties.VariableNames))
    fig0 = figure('Name','Frobenius norms','Color','w'); hold on; grid on;

    if any(T.method=="GRP")
        selA = (T.method=="GRP");
    else
        selA = (T.method==methods(1));
    end
    xA = double(T.batch(selA));
    yA = double(T.froX(selA));
    plot(xA, yA, '-', 'LineWidth',1.5, 'DisplayName','||A||_F');
    for i=1:numel(methods)
        sel = (T.method==methods(i)) & isfinite(T.froXhat);
        if any(sel)
            [bSorted, ord] = sort(double(T.batch(sel)));
            yhat = double(T.froXhat(sel));
            plot(bSorted, yhat(ord), '--', 'LineWidth',1.5, ...
                 'DisplayName',"||Â||_F ("+methods(i)+")");
        end
    end

    xlabel('Batch'); ylabel('Frobenius norm');
    title("Changes of norms over batches — " + titleSuffix);
    legend('Location','best'); hold off;

    saveas(fig0, fullfile(figDir, 'frobenius_norms_vs_batch.fig'));
    print(fig0, fullfile(figDir, 'frobenius_norms_vs_batch.png'), '-dpng','-r150');
else
    warning('Columns froX/froXhat not found in metrics table. Cannot plot Frobenius norms.');
end

fig1 = figure('Name','Reconstruction Error','Color','w');
hold on;
for i=1:numel(methods)
    sel = (T.method==methods(i)) & isfinite(T.recErr);
    if any(sel)
        [bSorted, ord] = sort(Tbatch(sel));
        y = double(T.recErr(sel));
        plot(bSorted, y(ord), '-o', 'LineWidth',1.5, 'MarkerSize',4, 'DisplayName',methods(i));
    end
end
grid on; xlabel('Batch'); ylabel('Reconstruction Error (Frobenius)');
title("Reconstruction Error vs. Batch — " + titleSuffix);
legend('Location','best'); hold off;
saveas(fig1, fullfile(figDir, 'recErr_vs_batch.fig'));
print(fig1, fullfile(figDir, 'recErr_vs_batch.png'), '-dpng','-r150');

fig2 = figure('Name','EVR','Color','w');
hold on;
for i=1:numel(methods)
    sel = (T.method==methods(i)) & isfinite(T.evr);
    if any(sel)
        [bSorted, ord] = sort(Tbatch(sel));
        y = double(T.evr(sel));
        plot(bSorted, y(ord), '-s', 'LineWidth',1.5, 'MarkerSize',4, 'DisplayName',methods(i));
    end
end
grid on; xlabel('Batch'); ylabel('Approx. EVR');
title("Explained Variance vs. Batch — " + titleSuffix);
legend('Location','best'); hold off;
saveas(fig2, fullfile(figDir, 'evr_vs_batch.fig'));
print(fig2, fullfile(figDir, 'evr_vs_batch.png'), '-dpng','-r150');

fig3 = figure('Name','Runtime','Color','w');
hold on;
for i=1:numel(methods)
    sel = (T.method==methods(i)) & isfinite(T.time_s);
    if any(sel)
        [bSorted, ord] = sort(Tbatch(sel));
        y = double(T.time_s(sel));
        plot(bSorted, y(ord), '-^', 'LineWidth',1.5, 'MarkerSize',4, 'DisplayName',methods(i));
    end
end
grid on; xlabel('Batch'); ylabel('Time (s)');
title("Runtime vs. Batch — " + titleSuffix);
legend('Location','best'); hold off;
saveas(fig3, fullfile(figDir, 'runtime_vs_batch.fig'));
print(fig3, fullfile(figDir, 'runtime_vs_batch.png'), '-dpng','-r150');

muRec = zeros(numel(methods),1);
muEvr = zeros(numel(methods),1);
muTim = zeros(numel(methods),1);
for i=1:numel(methods)
    sel = T.method==methods(i);
    muRec(i) = mean(double(T.recErr(sel)), 'omitnan');
    muEvr(i) = mean(double(T.evr(sel)), 'omitnan');
    muTim(i) = mean(double(T.time_s(sel)), 'omitnan');
end

fig4 = figure('Name','Mean RecErr','Color','w');
bar(muRec); set(gca,'XTick',1:numel(methods),'XTickLabel',methods);
ylabel('Mean RecErr'); title("Overall Mean RecErr — " + titleSuffix);
grid on; saveas(fig4, fullfile(figDir, 'mean_recErr.fig'));
print(fig4, fullfile(figDir, 'mean_recErr.png'), '-dpng','-r150');

fig5 = figure('Name','Mean EVR','Color','w');
bar(muEvr); set(gca,'XTick',1:numel(methods),'XTickLabel',methods);
ylabel('Mean EVR'); title("Overall Mean EVR — " + titleSuffix);
grid on; saveas(fig5, fullfile(figDir, 'mean_evr.fig'));
print(fig5, fullfile(figDir, 'mean_evr.png'), '-dpng','-r150');

fig6 = figure('Name','Mean Runtime','Color','w');
bar(muTim); set(gca,'XTick',1:numel(methods),'XTickLabel',methods);
ylabel('Mean Time (s)'); title("Overall Mean Runtime — " + titleSuffix);
grid on; saveas(fig6, fullfile(figDir, 'mean_runtime.fig'));
print(fig6, fullfile(figDir, 'mean_runtime.png'), '-dpng','-r150');

disp("نمودارها در پوشه ذخیره شد: " + figDir);

%% ==================== مرحله ۷: کشف الگوهای پرتکرار (FPM) ======================
disp("================================================================");
disp("مرحله 7: کشف الگوهای پرتکرار روی X و نسخه‌های بازسازی‌شده");
res7 = stage7_patterns(projectRoot, timestamp, struct('minSupport',minSupport,'maxLen',maxLen,'kRec_fd',kRec_fd));
fprintf('Pattern CSVs:\n%s\n%s\n%s\n%s\n', res7.orig, res7.grp, res7.ipca, res7.fd);

%% =========== مرحله ۸: تزریق ناهنجاری و ارزیابی تشخیص ===========================
disp("================================================================");
disp("مرحله 8: تزریق ناهنجاری و ارزیابی");
rep8 = stage8_anomaly(projectRoot, timestamp, b_inject, struct('nAnom',nAnom,'patternSize',patternSize,'kRec_fd',kRec_fd));
fprintf('Anomaly scores saved:\n%s\n%s\n%s\n', rep8.scores_grp, rep8.scores_ipca, rep8.scores_fd);

anomDir = fullfile(reportsDir,"anomaly");
injMat  = fullfile(anomDir, sprintf('batch_%03d_injected_%s.mat', b_inject, timestamp));
S8      = load(injMat, 'injIdx', 'Xb_mod');
injIdx  = S8.injIdx;
m_mod   = size(S8.Xb_mod,1);

scoreFiles = { ...
  fullfile(anomDir, sprintf('anomaly_scores_grp_b%03d_%s.csv',  b_inject, timestamp)), ...
  fullfile(anomDir, sprintf('anomaly_scores_ipca_b%03d_%s.csv', b_inject, timestamp)), ...
  fullfile(anomDir, sprintf('anomaly_scores_fd_b%03d_%s.csv',   b_inject, timestamp))};
labels = ["GRP","IPCA","FD"];

for si = 1:numel(scoreFiles)
    Tsc = readtable(scoreFiles{si});
    Tsc = sortrows(Tsc, 'row', 'ascend');
    scores = double(Tsc.score);
    rows   = double(Tsc.row);

    figA = figure('Name',"Anomaly Scores — " + labels(si),'Color','w');
    hold on; grid on;
    scatter(rows, scores, 12, 'filled', 'MarkerFaceAlpha', 0.25);
    showInjected = true;
    if showInjected
        scatter(injIdx, scores(injIdx), 22, 'r', 'filled');
    end
    xlabel('Row index in injected batch');
    ylabel('Anomaly score (row-wise residual)');
    title("Anomaly scores on Batch " + b_inject + " — " + labels(si) + " — " + datestr(now,'yyyy-mm-dd'));
    legend(["normal","injected"], 'Location','best');
    hold off;

    outFig = fullfile(reportsDir,'figures', ...
        sprintf('stage8_scores_%s_b%03d_%s', lower(labels(si)), b_inject, timestamp));
    saveas(figA, outFig + ".fig");
    print(figA, outFig + ".png", '-dpng','-r150');
end

%% ======================= مرحله ۹: گزارش نهایی و رتبه‌بندی =======================
disp("================================================================");
disp("مرحله 9: گزارش نهایی و رتبه‌بندی روش‌ها");
finalT = stage9_report(projectRoot, timestamp, struct('topN_patterns',50,'b_inject',b_inject));
disp("================================================================");
end
