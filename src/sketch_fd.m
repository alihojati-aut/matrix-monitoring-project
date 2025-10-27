function [B, Xhat] = sketch_fd(X_batch, ell, kRec)

    if nargin < 3 || isempty(kRec), kRec = 0; end

    Xb = double(full(X_batch));
    [m_b, n] = size(Xb);
    ell_eff = min([ell, m_b, n]);
    if ell_eff <= 0
        error('ell must be positive and <= min(size(X_batch)).');
    end

    [U, S, V] = svd(Xb, 'econ');

    s = diag(S);
    delta = (numel(s) >= ell_eff) * s(ell_eff)^2;
    S2 = max(S.^2 - delta, 0);
    S_shrink = sqrt(S2);

    C_shrunk = U * S_shrink * V.';
    if size(C_shrunk,1) >= ell
        B = C_shrunk(1:ell, :);
    else
        B = [C_shrunk; zeros(ell - size(C_shrunk,1), n)];
    end
    if kRec > 0
        [~, ~, Vb] = svd(B, 'econ');
        k = min([kRec, size(Vb,2)]);
        Vk = Vb(:, 1:k);
        Xhat = Xb * (Vk * Vk.');
    else
        Xhat = [];
    end
end
