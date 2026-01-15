% ichol with guarantee that H is spd

function L = icholspd(H, eta)
    if nargin < 2
        eta = 0.05;
    end

    opts.issym = true;
    opts.isreal = true;

    tic;
    lambda_min = eigs(H, 1, 'sm');
    minEig = toc
    tic;
    lambda_max = eigs(H, 1, 'lm');
    maxEig = toc

    droptol = eta * (lambda_min / lambda_max);

    ic_opts.type    = 'ict';
    ic_opts.droptol = droptol;
    ic_opts.michol  = 'off';

    tic;
    L = ichol(H, ic_opts);
    ichol=toc
end
