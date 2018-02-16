function [control,pic,T,p,constrained,squared,rho,lambda,verbose,phi0,u0,...
    accuracy,eps,tol,max_iter,max_cpu,online,warm_start,l2_prox,...
    if_trace_filter,last_part] = parse_input(y1,y2,control)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if iscolumn(y1) && iscolumn(y2)
    pic = 0;
elseif ismatrix(y1) && ismatrix(y2) ...
       && prod(size(y1) == size(y1')) && prod(size(y2) == size(y2'))
    pic = 1;
else
    error('y1 and y2 must be either column vectors or square matrices');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = length(y1)-length(y2); % y1 is the 'bandwidth signal' of length T+L+1
L = length(y2)-1; % y2 is the 'validation signal' of length L+1
% given y of length T+L+1, one can construct y1 = y and y2 = y(T+1:T+L+1)
if T < 0
    error('y1 must not be shorter than y2');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(control,'phi0') % both for MP and Nesterov
    control.phi0 = zeros(T+1,(T+1)^pic);
end
phi0 = control.phi0; % will be projected on l1-ball in the constrained case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(control,'u0') % only for MP
    control.u0 = zeros(L+1,(L+1)^pic);
end
u0 = control.u0; % will be projected onto l_q-ball
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = control.p;
verbose = control.verbose;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(control,'constrained')
    control.constrained = 1;
end
constrained = control.constrained;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if constrained % rho needed
    control.lambda = [];
    if isfield(control,'rho')
        if ~(control.rho >= 0)
            error('rho must be a non-negative number');
        end
    else
        error('constrained fitting, rho (filter norm) must be given');
    end
else % lambda needed
    control.rho = [];
    if isfield(control,'lambda')
        if ~(control.lambda > 0)
            error('lambda must be a strictly positive number');
        end
    else
        error('regularized fitting, lambda must be given');
    end
end
%% only for cvx
if ~isfield(control,'squared')
    control.squared = 0;
end
if control.squared && ~(~constrained && p == 2)
    error('min(|*|_p^2) only in the unconstrained case with p = 2');
end
squared = control.squared;  
%%
rho = control.rho;
lambda = control.lambda;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(control,'accuracy') % accuracy mode: 'abs' or 'rel'
    control.accuracy = 'abs';
end
accuracy = control.accuracy; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(control,'eps') % precision, absolute or relative
    error('eps (precision) not given');
end
eps = control.eps; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(control,'tol') % numerical tolerance for comparisons
    control.tol = 1e-8;
end
tol = control.tol;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(control,'max_iter') % maximum number of iterations
    control.max_iter = 100000;
end
max_iter = control.max_iter;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(control,'max_cpu') % maximum Mirror Prox computation time
    control.max_cpu = 100000;
end
max_cpu  = control.max_cpu;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(control,'online') % use stepsize search heuristics if true
    control.online = 1;
end
online = control.online;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(control,'warm_start') % initialize by signal DFT subspace
    control.warm_start = 0;
end
warm_start = control.warm_start;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sample = control.sample; % how many points to sample for approx Lip const
% memory = control.memory; % exponential multiplier for moving average
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(control,'l2_prox')
    control.l2_prox = 1;
end
l2_prox = control.l2_prox; % use l2 prox everywhere
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(control,'last_part')
    control.last_part = 50;
end
last_part = control.last_part; % heuristic, average over the last iters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(control,'if_trace_filter')
    control.if_trace_filter = 0;
end
if_trace_filter = control.if_trace_filter; % filters computed after each it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if verbose
    control
end