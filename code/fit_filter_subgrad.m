function sol = fit_filter_subgrad(y1,y2,control)
[control,pic,T,filt0,p,constrained,rho,lambda,verbose,accuracy,eps,tol,...
    max_iter,max_cpu,online,warm_start,l2_prox]...
        = parse_input(y1,y2,control);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Subgradient Descent %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b = dft(y2);
L = length(y2)-1;
% set prox parameters
[r, gamma] = prox_setup(l2_prox,(T+1)^(pic+1));
% currently only absolute accuracy, constrained problem
if p ~= inf
    error('Projected Subgradient Descent implemented only for p = inf');
end
time = cputime; % set timer for the Lipschitz constant calculation
% preprocessing : computing the Lipschitz constant
if verbose
    fprintf('\tConstructing the operator...\n');
end
A = construct_operator(T, rho, y1);
if l2_prox
    Lip = max(sqrt(sum(abs(A).^2,2))); % |A|_{2,inf} = max row l2-norm
else
    Lip = max(max(abs(A))); % |A|_{1,inf}
end
if verbose
    fprintf('\tDone. Time spent: %2.2f.\n', cputime-time);
end
time = cputime; % set timer for iterations
% maximum number of iters to do guaranteed by the theory
itersToDo = (Lip/eps)^2;
if online
    error('Adaptive stepsize not implemented for PSD');
else
    % non-adaptive diverging steps as recommended by theory
    step = 1/Lip * ones(max_iter,1) ./ sqrt(1:max_iter)';
    %step = 1/Lip * ones(max_iter,1) / sqrt(itersToDo);
end
% show info
if verbose
    fprintf('\nCalling the Subgradient Descent\n');
    fprintf('---------------------------------------\n');
    fprintf('Problem\n');
    fprintf('\tT, L\t\t:    %d, %d\n', T, L);
    fprintf('\teps\t\t:    %f\n',eps);
    fprintf('\tLipschitz const\t:    %f\n', Lip);
    fprintf('\tTheory # iters\t:    %d\n', itersToDo);
end
if ~warm_start
    phi = zeros(T+1,(T+1)^pic);
else
    phi = init_warm_start(pic,y1,b,T,L,rho,p);
end
phi_min = phi;
res = direct_operator(phi,rho,y1)-b;
f_min = norm(res(:),p);
iter = 0;
Pobj = f_min;
Cpu = cputime-time;
ifTerminate = (f_min < eps);
if verbose
    fprintf('\nOptimizer started.\n\n');
    fprintf(...
        strcat('ITER\t    OBJ\t\t    L1NORM\t    STEP\t    TIME\n\n'));
    fprintf('%d\t    ', 0);
    fprintf('%2.4f\t    ', f_min);
    fprintf('%2.4f\t    ', norm(phi_min(:),1));
    fprintf('%f\t    ', 0);
    fprintf('%2.2f\n', cputime-time);
end
%% Iterations
while ~ifTerminate
    %% Subgradient step
    res = direct_operator(phi, rho, y1) - b;
    [max_val,idx] = max(abs(res(:)));
    [i,j] = ind2sub(size(res),idx);
    u = zeros(size(b));
    u(i,j) = 1;
    g = sign(res(i,j)) * adjoint_operator(u,rho,y1);
    % subgradient step
    phi_plus = prox_grad(constrained,phi(:),r,gamma) - step(iter+1) * g(:);
    % must be constrained = 1
    phi = l1_lr_prox(constrained,phi_plus,r,gamma,tol);
    if pic, phi = vec2mat(phi); end
    % update minimal value
    res = direct_operator(phi,rho,y1)-b;
    fphi = norm(res(:),inf);
    %%
    if fphi < f_min
        phi_min = phi;
        f_min = fphi;
    end
    if verbose && ~mod(iter+1,verbose)
        fprintf('%d\t    ', iter+1);
        fprintf('%2.4f\t    ', f_min);
        fprintf('%2.4f\t    ', norm(phi_min(:),1));
        fprintf('%f\t    ', step(iter+1));
        fprintf('%2.2f\n', cputime-time);
    end
    Pobj = [Pobj;f_min];
    cpu = cputime-time;
    Cpu = [Cpu;cpu];
    ifCpu = (cpu >= max_cpu);
    ifIter = (iter+1 >= max_iter);
    ifTerminate = ifCpu || ifIter || (iter + 1 >= itersToDo);
    iter = iter+1;
end
%% Output
sol.status = 0;
if ifIter, sol.status = 1; end
if ifCpu, sol.status = 2; end
sol.cpu = cputime-time;
sol.iter = iter;
sol.online_iter_rate = [];
sol.online_stepsize_ratio = [];
sol.phi = phi_min;
sol.u = [];
sol.pobj = f_min;
sol.dobj = [];
sol.gap = [];
sol.filter = rho / sqrt(T+1)^(pic+1) * idft(phi_min);
% stats
sol.Pobj = Pobj;
sol.Dobj = [];
sol.Gap = [];
sol.Cpu = Cpu;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end