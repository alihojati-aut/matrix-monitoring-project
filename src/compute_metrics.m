function M = compute_metrics(X_batch, varargin)

    X = double(full(X_batch));
    froX2 = sum(X(:).^2);
    froX  = sqrt(froX2);

    p = inputParser;
    addParameter(p, 'recon',   []);
    addParameter(p, 'sketchR', []);
    addParameter(p, 'eps',     eps);
    parse(p, varargin{:});
    eps0 = p.Results.eps;

    if ~isempty(p.Results.recon)
        Xhat = double(p.Results.recon);
    elseif ~isempty(p.Results.sketchR)
        cellin = p.Results.sketchR;
        S = double(cellin{1}); R = double(cellin{2});
        Xhat = S * R.';
    else
        error('Provide either ''recon'', or ''sketchR'' (e.g. {S,R}).');
    end

    if ~isequal(size(Xhat), size(X))
        error('compute_metrics: size(Xhat) ~= size(X).');
    end

    Diff    = X - Xhat;
    diffF2  = sum(Diff(:).^2);
    diffF   = sqrt(diffF2);
    froXhat2 = sum(Xhat(:).^2);

    mu  = mean(X, 1);
    Xc  = X - mu;
    den = sum(Xc(:).^2);
    if den <= eps0
        den = max(froX2, eps0);
    end
    evr = 1 - (diffF2 / den);
    evr = max(0, min(1, evr));

    M = struct();
    M.froX     = froX;
    M.froX2    = froX2;
    M.froXhat  = sqrt(froXhat2);
    M.recErr   = diffF / max(froX, eps0);
    M.evr      = evr;
end
