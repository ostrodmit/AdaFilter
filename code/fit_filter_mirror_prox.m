function sol = fit_filter_mirror_prox(y1,y2,control)
[control,pic,T,p,constrained,squared,rho,lambda,verbose,phi0,u0,...
    accuracy,eps,tol,max_iter,max_cpu,online,warm_start,l2_prox,...
    if_trace_filter,last_part] = parse_input(y1,y2,control);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% Mirror Prox %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~constrained
        rho = 1; % currently unconstrained algo does not work properly with online = 1, especially with l2_prox = 1 -- probably due to numerical problems
    end
    b = dft(y2);
    L = length(y2)-1;
    if verbose
        fprintf('\nCalling Mirror Prox\n');
        fprintf('-----------------------------------------------\n');
        fprintf('Problem\n');
        fprintf('\tT, L\t\t:   %d, %d\n', T, L);
        %fprintf('\tDimension\t:   %d\n', 2*(T^(pic+1)+L^(pic+1)));
    end
    time = cputime; % set timer for the Lipschitz constant calculation
    if verbose
        disp('Computing the Lipschitz constant');
    end
    if l2_prox
        if verbose, fprintf('\tPower iteration...\n'); end
        nrep = 10;
        Lip = operator_norm_2_2(pic,y1,T,rho,nrep);
    else
        % use the diagonalisation upper bound
        if verbose, fprintf('\tUse the "almost diagonal" bound...\n'); end
        Lip = lip_const_upper_bound_1_p(y1,T,L,rho,pic);
%         % explicitly construct A and compute the Lipschitz constant
%         if verbose, fprintf('\tConstructing the operator...\n'); end
%         A = construct_operator(T, rho, y1);
%         if p == inf
%             Lip = max(max(abs(A))); %|A|_{1,inf}
%         else
%             Lip = max(sqrt(sum(abs(A).^2,1))); %|A|_{1,2} = max col l2-norm
%         end
    end
    theoryStep = 1/Lip;
    step = theoryStep * ones(max_iter, 1);
%     if strcmp(accuracy, 'abs')
%         maxItersToDo = round(exp(1)*(log((L+1)^(pic+1)+1) ...
%           + log((T+1)^(pic+1)+1)) / eps / step(1)) + 1;
%     else
%         maxItersToDo = round(exp(1)*(log((L+1)^(pic+1)+1) ...
%           + log((T+1)^(pic+1)+1)) / (eps * norm(b(:), p))/step(1)) + 1;
%     end
    if verbose
        fprintf('\tInverse Lip\t:   %0.4f\n', step(1));
%         if constrained
%             fprintf('\tTheory max # of iters\t:   %d\n', maxItersToDo);
%         end
    end
    
    % initial values for the last part heuristic
    m = last_part;
    if last_part
        last_part_frac = (1/m) : (1/m) : 1;
    end
    ss = zeros(m+1,1);
    s = zeros(m+1,1);
    last_part_pobj = zeros(m+1,1);
    last_part_dobj = zeros(m+1,1);
    last_part_gap = zeros(m+1,1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % phi
    if ~warm_start
        phi = phi0;
        if constrained && norm(phi(:),1) > 1 % normally <= 1
            phi = phi / norm(phi(:),1);
        end
    else % start from the signal DFT
        phi = init_warm_start(pic,y1,b,T,L,rho,p);
%        % find dual variable
%         if p == 2
%             u = conj(res) ./ norm(res(:), 2);
%         else % p == inf
%             [umax idx] = max(abs(res(:)));
%             u = zeros((L+1)^(pic+1),1);
%             u(idx) = conj(res(idx))/umax;
%             if pic, u = vec2mat(u); end
%         end
%        adj_u = adjoint_operator(u, rho, y1);
%        dobj = max( 0, - norm(adj_u, inf)... 
%                       - real(dot(u(:),b(:))) );
    end
    % u
    u = u0;
    if p == 2, q = 2; else q = 1; end
    if norm(u(:),q) > 1 % normally <= 1
        u = u / norm(u(:),q);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nabla_phi = zeros(T+1,(T+1)^pic);    
    mean_nabla_pphi = zeros(T+1,(T+1)^pic, m+1);
    mean_pphi = zeros(T+1,(T+1)^pic, m+1);
    
    nabla_u = zeros(L+1,(L+1)^pic);
    mean_nabla_uu = zeros(L+1,(L+1)^pic, m+1);
    mean_uu = zeros(L+1,(L+1)^pic, m+1);
    
    last_part_mean_nabla_pphi = mean_nabla_pphi;
    last_part_mean_pphi = mean_pphi;
    last_part_mean_nabla_uu = mean_nabla_uu;
    last_part_mean_uu = mean_uu;
    
    %%%%%%%%%%%%%%%%%%%% Intializing parameters for DGF %%%%%%%%%%%%%%%%%%%
    % set prox parameters
    [r_phi,gamma_phi] = prox_setup(l2_prox,(T+1)^(pic+1));
    [r_u,gamma_u] = prox_setup(l2_prox||p==2,(L+1)^(pic+1));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    iter = 0;
    best_k = 0;
    next_iter = ones(m+1,1);
    int_iter = 0;
    %old_int_iter = 0;
    res = direct_operator(phi, rho, y1) - b;
    pobj = norm(res(:),p);
    dobj = 0;
    gap = pobj-dobj;
    relErrUb = gap/dobj;
    cpu = cputime - time;
    % collect stats
    Pobj = pobj;
    Dobj = dobj;
    Gap = gap;
    RelErrUb = relErrUb;
    Cpu = cpu;
    Int_iter = int_iter;
    % trace computed filter
    if if_trace_filter && ~pic,
        new_filter = squeeze(rho / sqrt(T+1) * idft(phi));
        Filter = new_filter;
    end
    % check terminal condition
    ifCpu = (cpu >= max_cpu);
    ifIter = (iter >= max_iter);
    ifTerminate = ifCpu || ifIter ...
        || (strcmp(accuracy,'abs') && gap < eps) ...
        || (strcmp(accuracy,'rel') && relErrUb < eps);
    status = 0;
    if ifCpu, status = 1; end
    if ifIter, status = 2; end
    if verbose
        fprintf('\nOptimizer started.\n\n');
        fprintf(strcat('ITER\t     INT\t    POBJ\t    DOBJ\t    GAP\t\t      ', ...
                       '    REL\t\t    STEP\t    TIME\n\n'));
        fprintf('%4d\t    ', iter);
        fprintf('%4d\t    ', int_iter);
        fprintf('%2.4f\t    ', pobj);
        fprintf('%2.4f\t    ', dobj);
        fprintf('%2.4f\t    ', gap);
        fprintf('%0.2f\t    ', relErrUb);
        fprintf('%1.2f\t    ', 0);
        fprintf('%2.2f\n', cputime-time);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Main loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while ~ifTerminate
        %%%%%%%%%%%%%%%%%%%%%%%%%%% make step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~online || iter <= 1 % use the fixed stepsize
            [phi, u, pphi, uu, nabla_pphi, nabla_uu]...
                = take_step_mirror_prox(constrained,phi,u,y1,b,rho,... 
                [r_phi,r_u],[gamma_phi,gamma_u],pic,step(iter+1),p,tol,lambda);
            %old_int_iter = int_iter;
            int_iter = int_iter + 1;
        else % adaptive stepsize heuristic
            q_inc = 2;
            q_dec = sqrt(3);
            if step(iter) >= step(iter-1)
                step(iter+1) = step(iter) * q_inc;
            else
                step(iter+1) = step(iter);
            end
            % save current state
            phi_init = phi;
            u_init = u;
            nabla_phi_init = nabla_phi;
            nabla_u_init = nabla_u;
            % step for the first time
            [phi, u, pphi, uu, nabla_pphi, nabla_uu,...
                prox_grad_phi_ext, prox_grad_phi_old]...
                = take_step_mirror_prox(constrained,phi_init,u_init,y1,b,...
                rho,[r_phi,r_u],[gamma_phi,gamma_u],pic,step(iter+1),p,tol,lambda);
            %old_int_iter = int_iter;
            int_iter = int_iter + 1;
            while ~mirror_prox_stepsize_invariant(constrained, ...
                    phi_init(:), pphi(:), phi(:), u_init(:), uu(:), u(:),...
                    nabla_phi_init(:), nabla_pphi(:), nabla_uu(:),...
                    prox_grad_phi_old(:), prox_grad_phi_ext(:),...
                    [r_phi r_u], [gamma_phi gamma_u], step(iter+1),...
                    theoryStep,tol)
                step(iter+1) = step(iter+1)/q_dec;
                [phi, u, pphi, uu, nabla_pphi, nabla_uu,...
                    prox_grad_phi_ext]...
                    = take_step_mirror_prox(constrained,phi_init,u_init,y1,b,...
                    rho,[r_phi,r_u],[gamma_phi,gamma_u],pic,step(iter+1),p,tol,lambda);
                int_iter = int_iter + 1;
            end
        end
        %%%%%%%%%%%%%%%%%% prepare for the next iteration %%%%%%%%%%%%%%%%%
        st = step(iter+1);
        % save old sum of weights
        ss(1) = s(1);
        % update sum of weights
        s(1)  = sum(step(1:iter+1));
        % update averages
        mean_nabla_pphi(:,:,1) = (mean_nabla_pphi(:,:,1) * ss(1) ...
            + nabla_pphi * st) / s(1);
        mean_nabla_uu(:,:,1) = (mean_nabla_uu(:,:,1) * ss(1) ...
            + nabla_uu * st) / s(1);
        mean_pphi(:,:,1) = (mean_pphi(:,:,1) * ss(1) ...
            + pphi * st) / s(1);
        mean_uu(:,:,1) = (mean_uu(:,:,1) * ss(1) ...
            + uu * st) / s(1);
        last_part_mean_nabla_pphi(:,:,1) = mean_nabla_pphi(:,:,1);
        last_part_mean_pphi(:,:,1) = mean_pphi(:,:,1);
        last_part_mean_nabla_uu(:,:,1) = mean_nabla_uu(:,:,1);
        last_part_mean_uu(:,:,1) = mean_uu(:,:,1);
        %%%%%%%%%%%%%% re-evaluate primal and dual objective %%%%%%%%%%%%%%
        l_p_m_n_uu = last_part_mean_nabla_uu(:,:,1); 
        l_p_m_n_pphi  = last_part_mean_nabla_pphi(:,:,1);
        l_p_m_uu = last_part_mean_uu(:,:,1);
        last_part_pobj(1) = norm(l_p_m_n_uu(:),p);
%         last_part_dobj(1) = -norm(l_p_m_n_pphi(:),inf) ...
%             - real(dot(l_p_m_uu(:),b(:)));
        last_part_dobj(1) = max(0,-norm(l_p_m_n_pphi(:),inf) ...
            - real(dot(l_p_m_uu(:),b(:)))); % can only be positive        
        if ~constrained
            l_p_m_pphi = last_part_mean_pphi(:,:,1);
            last_part_pobj(1) = last_part_pobj(1) + lambda * norm(l_p_m_pphi(:),1);
            last_part_dobj(1) = 0;
        end
        %%%%%%%%%%%%%%%%%%%% re-evaluate the dual gap %%%%%%%%%%%%%%%%%%%%%
        last_part_gap(1) = ... 
            last_part_pobj(1) - last_part_dobj(1);
        pobj = last_part_pobj(1); dobj = last_part_dobj(1);
        gap = last_part_gap(1);
        % other divisors
        best_k = 0;
        for k = 1:m
            if iter == next_iter(k+1)
                %%%%%%%%%%%%%% re-calculate running averages %%%%%%%%%%%%%%
                % save old sum of weights
                ss(k+1) = s(k+1);
                % update sum of weights
                s(k+1)  = s(1);
                % update last part averages
                last_part_mean_nabla_pphi(:,:,k+1) = ...
                   ( mean_nabla_pphi(:,:,1) * s(k+1) ... 
                    -mean_nabla_pphi(:,:,k+1) * ss(k+1))/(s(k+1)-ss(k+1));
                last_part_mean_nabla_uu(:,:,k+1) = ...
                   ( mean_nabla_uu(:,:,1) * s(k+1) ... 
                    -mean_nabla_uu(:,:,k+1) * ss(k+1))/(s(k+1)-ss(k+1));
                last_part_mean_pphi(:,:,k+1) = ...
                   ( mean_pphi(:,:,1) * s(k+1) ... 
                    -mean_pphi(:,:,k+1) * ss(k+1))/(s(k+1)-ss(k+1));
                last_part_mean_uu(:,:,k+1) = ...
                   ( mean_uu(:,:,1) * s(k+1) ... 
                    -mean_uu(:,:,k+1) * ss(k+1))/(s(k+1)-ss(k+1));
                % update full averages
                mean_nabla_pphi(:,:,k+1) = mean_nabla_pphi(:,:,1);
                mean_nabla_uu(:,:,k+1) = mean_nabla_uu(:,:,1);
                mean_pphi(:,:,k+1) = mean_pphi(:,:,1);
                mean_uu(:,:,k+1) = mean_uu(:,:,1);
                %%%%%%%%%% re-evaluate primal and dual objective %%%%%%%%%%
                l_p_m_n_uu = last_part_mean_nabla_uu(:,:,k+1); 
                l_p_m_n_pphi  = last_part_mean_nabla_pphi(:,:,k+1);
                l_p_m_uu = last_part_mean_uu(:,:,k+1);
                last_part_pobj(k+1) = norm(l_p_m_n_uu(:),p);
%                 last_part_dobj(k+1) = - norm(l_p_m_n_pphi(:),inf) ... 
%                        - real(dot(l_p_m_uu(:),b(:)));
                last_part_dobj(k+1) = max(0,-norm(l_p_m_n_pphi(:),inf) ... 
                       - real(dot(l_p_m_uu(:),b(:)))); % can only be positive
                if ~constrained
                    l_p_m_pphi = last_part_mean_pphi(:,:,k+1);
                    last_part_pobj(k+1) = last_part_pobj(k+1)... 
                                          + lambda * norm(l_p_m_pphi(:),1);
                    last_part_dobj(k+1) = 0;
                end
                %%%%%%%%%%%%%%%% re-evaluate the dual gap %%%%%%%%%%%%%%%%%
                last_part_gap(k+1) = ... 
                    last_part_pobj(k+1) - last_part_dobj(k+1);
                if last_part_gap(k+1) < gap
                    best_k = k;
                    pobj = last_part_pobj(k+1); dobj = last_part_dobj(k+1);
                    gap = last_part_gap(k+1);
                end
                next_iter(k+1) = ceil(iter / (1-last_part_frac(k)));
            end
        end
        % take into account the last iteration (only in the dual)
        % dobj = max(dobj, -norm(nabla_pphi(:),inf)+real(dot(nabla_uu(:),b(:))));
        if ~constrained
            dobj = 0; % no precision control in the unconstrained case
        end
        relErrUb = gap/dobj;
        cpu = cputime - time;
        % Collect stats
        Pobj = [Pobj;pobj];
        Dobj = [Dobj;dobj];
        Gap = [Gap;gap];
        RelErrUb = [RelErrUb;relErrUb];
        Cpu = [Cpu;cpu];
        Int_iter = [Int_iter;int_iter];
        if if_trace_filter && ~pic,
            new_filter ... 
                = squeeze(rho / sqrt(T+1) * idft(mean_pphi(:,:,best_k+1)));
            Filter = [Filter new_filter];
        end
        % check terminal condition
        ifCpu = (cpu >= max_cpu);
        ifIter = (iter+1 >= max_iter);
        ifTerminate = ifCpu || ifIter...
            || strcmp(accuracy,'abs') && gap < eps...
            || strcmp(accuracy,'rel') && relErrUb < eps;...
            %|| iter > maxItersToDo && constrained;
        status = 0;
        if ifIter, status = 1; end
        if ifCpu, status = 2; end
        if verbose && ~mod(iter+1,verbose)
            fprintf('%4d\t    ', iter+1);
            fprintf('%4d\t    ', int_iter);
            fprintf('%2.4f\t    ', pobj);
            fprintf('%2.4f\t    ', dobj);
            fprintf('%2.4f\t    ', gap);
            fprintf('%0.2f\t    ', relErrUb);
            fprintf('%1.2f\t    ', st);
            fprintf('%2.2f\n', cputime-time);
        end
        iter = iter + 1;
    end
    if verbose
        fprintf('\nOptimizer terminated.  Time: %0.2f\n\n', cputime-time);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sol.status = status; % 1 if max_iter, 2 if max_cpu reached, 0 otherwise
    sol.cpu = cputime-time; % total running time in sec
    sol.iter = iter; % total number of Mirror Prox iterations done
    % stepsize heuristic efficiency : ratio of iteration counts
    sol.online_iter_rate = int_iter/iter;
    % stepsize heuristic efficiency : theory to practice stepsize ratio
    sol.online_stepsize_ratio = mean(step(1:iter))/theoryStep;
    if iter > 0
        sol.phi = mean_pphi(:,:,best_k+1); % reported primal solution
        sol.u = mean_uu(:,:,best_k+1); % reported dual solution        
    else
        sol.phi = phi;
        sol.u = u;
    end
    sol.pobj = pobj; % primal objective value at the reported solution
    sol.dobj = dobj; % dual objective value at the reported solution
    sol.lowerb = []; % in Nesterov it is instead of dobj
    sol.gap = gap; % primal-dual gap
    % filter
    sol.filter = rho / sqrt(T+1)^(pic+1) * idft(sol.phi);
    % stats
    sol.Pobj = Pobj;
    sol.Dobj = Dobj;
    sol.Lowerb = [];
    sol.RelErrUb = RelErrUb;
    sol.Gap = Gap;
    sol.Cpu = Cpu;
    sol.Int_iter = Int_iter;
    if if_trace_filter && ~pic,
        sol.Filter = Filter;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

% function ifTrue = mirror_prox_stepsize_invariant(constrained,...%old ver
%     phi_old, phi_ext, phi_new, u_old, u_ext, u_new, nabla_phi_old,...
%     nabla_phi_ext, nabla_u_old, nabla_u_ext, r, gamma, st, theoryStep,...
%     tol)
% r_phi = r(1); r_u = r(2); gamma_phi = gamma(1); gamma_u = gamma(2);
% D_new_ext = bregman_div(constrained, phi_new, phi_ext, r_phi, gamma_phi)...
%     + bregman_div(1, u_new, u_ext, r_u, gamma_u);
% D_ext_old = bregman_div(constrained, phi_ext, phi_old, r_phi, gamma_phi)...
%     + bregman_div(1, u_ext, u_old, r_u, gamma_u);
% z_ext = [phi_ext; u_ext];
% z_new = [phi_new; u_new];
% g_ext = [nabla_phi_ext; -nabla_u_ext];
% g_old = [nabla_phi_old; -nabla_u_old];
% delta = st * real(dot(g_ext-g_old, z_ext-z_new)) - D_new_ext - D_ext_old;
% ifTrue = (delta <= tol);% ifTrue = (delta <= 0) || (st <= theoryStep);
% end

function ifTrue = mirror_prox_stepsize_invariant(constrained,...%new ver
    phi_old, phi_ext, phi_new, u_old, u_ext, u_new, nabla_phi_old,...
    nabla_phi_ext, nabla_u_ext, prox_grad_phi_old, prox_grad_phi_ext,...
    r, gamma, st, theoryStep, tol)
% condition to check in the adaptive stepsize mode
r_phi = r(1); r_u = r(2); gamma_phi = gamma(1); gamma_u = gamma(2);
% maxel = max([abs(phi_new); abs(phi_old); abs(phi_ext);... 
%              abs(u_new); abs(u_old); abs(u_ext);... 
%              abs(nabla_phi_ext); abs(nabla_u_ext)]);
% phi_new_temp = phi_new / maxel; 
% phi_old_temp = phi_old / maxel;
% D_phi = bregman_div(constrained, phi_new_temp, phi_old_temp, r_phi, gamma_phi);
D_phi = bregman_div(constrained, phi_new, phi_old, r_phi, gamma_phi);
% u_new_temp = u_new / maxel;
% u_old_temp = u_old / maxel;
% D_u = bregman_div(1, u_new_temp, u_old_temp, r_u, gamma_u);
D_u = bregman_div(1, u_new, u_old, r_u, gamma_u);
D_new_old = D_phi + D_u;
% D_new_old = D_new_old * maxel^2;
z_ext = [phi_ext; u_ext];
z_new = [phi_new; u_new];
g = st * [nabla_phi_ext; -nabla_u_ext];
% z_ext_temp = z_ext / maxel;
% z_new_temp = z_new / maxel;
% g_ext_temp = g_ext / maxel;
% delta = st * maxel^2 * real(dot(g_ext_temp, z_ext_temp-z_new_temp)) - D_new_old;
delta = real(dot(g,z_ext-z_new))-D_new_old;
if ~constrained
    subg = prox_grad_phi_old - prox_grad_phi_ext - st * nabla_phi_old;
    delta = delta + real(dot(subg,phi_ext-phi_new));
end
ifTrue = (delta <= tol);
%ifTrue = (delta <= tol) || (st <= theoryStep);
end