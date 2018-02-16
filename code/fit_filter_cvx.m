function sol = fit_filter_cvx(y1,y2,control)
[control,pic,T,p,constrained,squared,rho,lambda,verbose]...
    = parse_input(y1,y2,control);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CVX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~constrained
        rho = 1;
    end
    b = dft(y2);
    time = cputime;
    if verbose
        fprintf('Constructing the operator...\n');
    end
    A = construct_operator(T, rho, y1);
    if verbose
        fprintf('Done. Time spent: %2.2f.\n', cputime-time);
    end
    time = cputime; % set timer for the iterations
    if verbose
        cvx_begin
    else
        cvx_begin quiet
    end
    % cvx_solver mosek
    cvx_precision best
        variable phi((T+1)^(pic+1)) complex
        if constrained
            minimize ( norm (b(:) - A * phi, p) )
            subject to
                norm(phi,1) <= 1;
        else
            if ~squared
                minimize (norm (b(:) - A * phi, p) + lambda * norm(phi,1))
            else % warning, quadratic form!
                minimize (square_pos(norm (b(:) - A * phi, p)) + lambda * norm(phi,1))
            end
        end
    cvx_end
    if pic
        phi = vec2mat(phi);
    end
    if strcmp(cvx_status, 'Solved')
        sol.status = 0;
    else
        sol.status = 1;
    end
    sol.cpu = cputime-time;
    sol.iter = [];
    sol.online_iter_rate = [];
    sol.online_stepsize_ratio = [];
    sol.phi = phi;
    sol.u = [];
    sol.pobj = cvx_optval;
    sol.dobj = [];
    sol.gap = [];
    sol.filter = rho / sqrt(T+1)^(pic+1) * idft(phi);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end