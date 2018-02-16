function sol = fit_filter_nesterov(y1,y2,control)
[control,pic,T,p,constrained,squared,rho,lambda,verbose,phi0,u0,...
    accuracy,eps,tol,max_iter,max_cpu,online,warm_start,l2_prox,...
    if_trace_filter,last_part] = parse_input(y1,y2,control);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Nesterov's Accelerated Gradient Descent %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~constrained
    rho = 1; % if unconstrained, penalize the squared norm!
end
b = dft(y2);
L = length(y2)-1;
% set prox parameters
[r, gamma] = prox_setup(l2_prox,(T+1)^(pic+1));
% currently only absolute accuracy (there is no gap), constrained problem
if p ~= 2
    error('The Nesterov method implemented only for p = 2');
end
time = cputime; % set timer for the Lipschitz constant calculation
if verbose
    disp('Computing the Lipschitz constant');
end
% preprocessing : computing the Lipschitz constant
if l2_prox %|2A^*A|_{2,2}=2|A|_{2,2}^2
    if verbose, fprintf('\tPower iteration...\n'); end
    nrep = 10;
    Lip = 2*(operator_norm_2_2(pic,y1,T,rho,nrep))^2;
else  % |2A^*A|_{1,inf}=2|A|_{1,2}^2
    % use the diagonalisation upper bound
    if verbose, fprintf('\tUse the "almost-diagonal" bound...\n'); end
    Lip = 2*(lip_const_upper_bound_1_p(y1,T,L,rho,pic))^2;
%     % explicitly construct A and compute the Lipschitz constant
%     if verbose, fprintf('\tConstructing the operator...\n'); end
%     A = construct_operator(T, rho, y1);
%     Lip = 2*max(sum(abs(A).^2,1));
end
if online
    LLip = Lip;
end
if verbose
    fprintf('\tDone. Time spent: %2.2f\n', cputime-time);
end
time = cputime; % set timer for iterations
% % maximum number of iters to do guaranteed by the theory
% if strcmp(accuracy,'abs')
%     maxItersToDo = 2*sqrt(Lip)/eps; % eps is the error for f_*, not f_*^2
% else
%     %not implemented
% end
% show info
if verbose
    fprintf('\nCalling the Nesterov method\n');
    fprintf('---------------------------------------\n');
    fprintf('Problem\n');
    fprintf('\tT, L\t\t:    %d, %d\n', T, L);
    fprintf('\teps\t\t:    %f\n',eps);
    fprintf('\tLipschitz const\t:    %f\n', Lip);
%     if constrained 
%         fprintf('\tTheory max # iters\t:    %d\n', maxItersToDo);
%     end
end
if ~warm_start
    phi = phi0;
    if constrained && norm(phi(:),1) > 1 % normally <= 1
        phi = phi / norm(phi(:),1);
    end
else % start from the signal DFT
    phi = init_warm_start(pic,y1,b,T,L,rho,p);
end
A = 0; % sum of weights
Gx = zeros(size(phi));
domega0 = prox_grad(constrained,phi(:),r,gamma); % weighted sum of gradients will be added
if pic, domega0 = vec2mat(domega0); end
phi_min = phi;
res = direct_operator(phi,rho,y1)-b;
fphi = norm(res(:),2);
if constrained
    obj = fphi;
else
    obj = fphi^2 + lambda * norm(phi_min(:),1);
end
best_obj = obj;
lowerb = 0;
gap = best_obj-lowerb;
Pobj = best_obj;
Lowerb = lowerb;
LinAggr = 0;
Gap = gap;
relErrUb = inf;
RelErrUb = relErrUb;
Cpu = cputime-time;
iter = 0;
int_iter = 0;
Int_iter = int_iter;
% trace computed filter
if if_trace_filter && ~pic,
    new_filter = squeeze(rho / sqrt(T+1) * idft(phi));
    Filter = new_filter;
end
if verbose
    fprintf('\nOptimizer started.\n\n');
    fprintf(...
        strcat('ITER\t     INT\t    BOBJ\t\t    LOW\t\t    GAP\t\t    REL\t\t    L1NORM\t    STEP\t    TIME\n\n'));
    fprintf('%4d\t    ', 0);
    fprintf('%4d\t    ', 0);
    fprintf('%2.4f\t    ', best_obj);
    fprintf('%2.4f\t    ', lowerb);
    fprintf('%2.4f\t    ', gap);
    fprintf('%2.4f\t    ', relErrUb);
    fprintf('%2.4f\t    ', norm(phi_min(:),1));
    fprintf('%f\t    ', 0);
    fprintf('%2.2f\n', cputime-time);
end
ifIter = 0;
ifCpu = 0;
ifTerminate = strcmp(accuracy,'abs') && gap < eps...
    || strcmp(accuracy,'rel') && relErrUb < eps;%...
%% Iterations
while ~ifTerminate
    LinAggr = LinAggr * A; % weighted sum of fx^2-<x,g>, for ave_lowerb
    %% Nesterov step
    if ~online
        [phi,x,x_hat,z,gx,gphi,a,A,fx,fphi]...
            = take_step_nesterov(pic,constrained,...
            y1,b,rho,Lip,phi,Gx-domega0,A,r,gamma,tol,lambda);
        int_iter = int_iter+1;
    else % online
        LLip = LLip / 2;
        % save current state
        phi_init = phi; A_init = A;
        % step for the first time
        [phi,x,x_hat,z,gx,gphi,a,A,fx,fphi]...
            = take_step_nesterov(pic,constrained,...
            y1,b,rho,LLip,phi_init,Gx-domega0,A_init,r,gamma,tol,lambda);
        int_iter = int_iter+1;
        while ~nesterov_step_invariant(phi(:),x(:),x_hat(:),z(:),...
                gx(:),fx,fphi,A,constrained,r,gamma,tol);
            LLip = LLip * sqrt(2);
            [phi,x,x_hat,z,gx,gphi,a,A,fx,fphi]...
                = take_step_nesterov(pic,constrained,y1,b,rho,...
                    LLip,phi_init,Gx-domega0,A_init,r,gamma,tol,lambda);
            int_iter = int_iter+1;
        end
    end
    Gx = Gx + a * gx;
    %%
    if constrained
        obj = fphi;
    else
        obj = fphi^2 + lambda * norm(phi(:),1);
    end
    % update the best objective value
    if obj < best_obj
        phi_min = phi;
        best_obj = obj;
    end
    %%
    if constrained
        LinAggr = (LinAggr + a * (fx^2 - real(dot(gx(:),x(:)))))/A;
        % the old value was already multiplied by A before the step
        ave_lowerb = LinAggr - norm(Gx(:)/A,inf);
        cur_lowerb_x = fx^2 - real(dot(gx(:),x(:))) - norm(gx(:),inf);
        cur_lowerb_phi = fphi^2 - real(dot(gphi(:),phi(:))) - norm(gphi(:),inf);
        cur_lowerb = max(cur_lowerb_x,cur_lowerb_phi);
        lowerb = max(lowerb,sqrt(max(0,max(ave_lowerb,cur_lowerb))));
    else
        lowerb = 0; % no precision control in the unconstrained case
    end
    gap = best_obj - lowerb;
    relErrUb = gap/lowerb;
    if verbose && ~mod(iter+1,verbose)
        fprintf('%4d\t    ', iter+1);
        fprintf('%4d\t    ', int_iter);
        fprintf('%2.4f\t    ', best_obj);
        fprintf('%2.4f\t    ', lowerb);
        fprintf('%2.4f\t    ', gap);
        fprintf('%2.4f\t    ', relErrUb);
        fprintf('%2.4f\t    ', norm(phi(:),1));
        fprintf('%f\t    ', a);
        fprintf('%2.2f\n', cputime-time);
    end
    Pobj = [Pobj;best_obj];
    Lowerb = [Lowerb;lowerb];
    Gap = [Gap;gap];
    RelErrUb = [RelErrUb;relErrUb];
    cpu = cputime-time;
    Cpu = [Cpu;cpu];
    Int_iter = [Int_iter;int_iter];
    if if_trace_filter && ~pic,
        new_filter ... 
            = squeeze(rho / sqrt(T+1) * idft(phi(:,:)));
        Filter = [Filter new_filter];
    end
    ifCpu = (cpu >= max_cpu);
    ifIter = (iter+1 >= max_iter);
    ifTerminate = ifCpu || ifIter...
        || strcmp(accuracy,'abs') && gap < eps...
        || strcmp(accuracy,'rel') && relErrUb < eps;%...
        %|| iter > maxItersToDo && constrained;
    iter = iter+1;
end
%% Output
sol.status = 0;
if ifIter, sol.status = 1; end
if ifCpu, sol.status = 2; end
sol.cpu = cputime-time;
sol.iter = iter;
sol.int_iter = int_iter;
sol.online_iter_rate = int_iter/iter;
sol.online_stepsize_ratio = [];
sol.phi = phi_min;
sol.u = [];
sol.pobj = best_obj;
sol.dobj = [];
sol.lowerb = lowerb;
sol.gap = gap;
sol.filter = rho / sqrt(T+1)^(pic+1) * idft(phi_min);
% stats
sol.Pobj = Pobj;
sol.Dobj = [];
sol.Lowerb = Lowerb;
sol.RelErrUb = RelErrUb;
sol.Gap = Gap;
sol.Cpu = Cpu;
sol.Int_iter = Int_iter;
if if_trace_filter && ~pic,
    sol.Filter = Filter;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function [phi,x,x_hat,z,gx,gphi,a,A,fx,fphi]...
    = take_step_nesterov(pic,constrained,y1,b,rho,Lip,phi,Gx,A,r,gamma,tol,lambda)
a = 1/2/Lip + sqrt(1/4/Lip^2+A/Lip);
A = A + a;
tau = a/A;
z = l1_lr_prox(constrained,-Gx(:),r,gamma,tol,A,lambda); % attention for the step!
if pic, z = vec2mat(z); end
x = tau * z + (1-tau) * phi;
xtemp = direct_operator(x, rho, y1) - b;
fx = norm(xtemp(:),2);
gx = 2 * adjoint_operator(xtemp, rho, y1);
%Gx = Gx + a * gx;
domegaz = prox_grad(constrained,z(:),r,gamma);
eta = domegaz - a*gx(:);
x_hat = l1_lr_prox(constrained,eta,r,gamma,tol,a,lambda);
if pic, x_hat = vec2mat(x_hat); end
phi = tau * x_hat + (1-tau) * phi;
phitemp = direct_operator(phi,rho,y1)-b;
fphi = norm(phitemp(:),2);
gphi = 2 * adjoint_operator(phitemp, rho, y1);
end

function ifTrue = nesterov_step_invariant(y,x,x_hat,z,gx,fx,fy,A,...
    constrained,r,gamma,tol)
delta = fx^2 - fy^2 + real(dot(gx,(y-x)))...
    + bregman_div(constrained,x_hat,z,r,gamma)/A;
ifTrue = delta >= -tol;
end