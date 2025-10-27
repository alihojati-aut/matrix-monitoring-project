function M = compute_metrics(X_batch, varargin)

    X = double(full(X_batch));
    froX2 = sum(X(:).^2);
    froX  = sqrt(froX2);

    p = inputParser;
    addParameter(p, 'recon', []);
    addParameter(p, 'sketchR', []);
    parse(p, varargin{:});

    if ~isempty(p.Results.recon)
        Xhat = double(p.Results.recon);
    elseif ~isempty(p.Results.sketchR)
        cellin = p.Results.sketchR;
        S = double(cellin{1}); R = double(cellin{2});
        Xhat = S * R.';
    else
        error('Provide either ''recon'', or ''sketchR'' (e.g. {S,R}).');
    end

    diffF = norm(X - Xhat, 'fro');
    froXhat2 = sum(Xhat(:).^2);

    M = struct();
    M.froX    = froX;
    M.froX2   = froX2;
    M.froXhat = sqrt(froXhat2);
    M.recErr  = diffF / max(froX, eps);
    M.evr     = min(froXhat2 / max(froX2, eps), 1);
end
