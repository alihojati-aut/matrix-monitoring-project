function [Z_batch, Xhat, ipca] = sketch_ipca(X_batch, kComponents, ipca)
    Xb = double(full(X_batch));

    mu = mean(Xb, 1);
    Xc = bsxfun(@minus, Xb, mu);

    [U,S,V] = svd(Xc, 'econ');
    k = min([kComponents, size(V,2)]);
    W = V(:, 1:k);

    Z_batch = Xc * W;
    Xhat    = Z_batch * W' + repmat(mu, size(Xb,1), 1);

    ipca = struct('Components', W', 'Mean', mu);
end
