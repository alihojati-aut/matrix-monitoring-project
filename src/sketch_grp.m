function [S_batch, Xhat, R] = sketch_grp(X_batch, dTarget, rngSeed)

    Xb = double(full(X_batch));

    nFeatures = size(Xb, 2);
    if nargin >= 3 && ~isempty(rngSeed), rng(rngSeed); end
    R = randn(nFeatures, dTarget) / sqrt(dTarget);

    S_batch = Xb * R;

    Xhat = S_batch * R.';
end
