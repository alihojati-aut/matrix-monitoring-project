function [S_batch, Xhat, R] = sketch_grp(X_batch, dTarget, rngSeed)
    X  = double(full(X_batch));
    [m,n] = size(X);

    mu = mean(X,1);
    Xc = bsxfun(@minus, X, mu);

    if nargin>=3 && ~isempty(rngSeed), rng(rngSeed); end
    Omega   = randn(n, dTarget) / sqrt(dTarget);
    S_batch = Xc * Omega;

    d = dTarget;

    lambda = 1e-3 * (norm(S_batch,'fro')^2 / max(m*d,1));
    Rt = (S_batch.' * S_batch + lambda*eye(d)) \ (S_batch.' * Xc);
    R  = Rt.';

    Xc_hat = S_batch * Rt;
    Xhat   = bsxfun(@plus, Xc_hat, mu);
end
