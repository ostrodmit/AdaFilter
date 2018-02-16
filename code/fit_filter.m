function sol = fit_filter(y1, y2, control)
%
% FIT_FILTER Fit denoising filter using two samples by convex programming.
%
% Input params
%   y1 : convolution sample (1d or 2d square) of length L+T+1
%   y2 : validation sample (1d or 2d square) of length L+1, i.e. no shorter
%        than y1.
%   control : structure of control parameters with the following fields
%           (we mark mandatory ones with *; in brackets we show the correct
%           value range and default value if any):
%       p [{2,inf},2] : l_2 or l_inf minimization. If constrained=0
%           and solver='nes', penalized l_2^2 objective will be minimized.
%       constrained [{0,1},1] : constrained or unconstrained problem.
%       *rho [{>0}] : filter norm constraint. Needed if contsrained = 1.
%       *lambda [{>0}] : regularization parameter. Needed if constrained =
%           0. Normalized on the *unit* l1 ball, i.e. before the change of
%           variables the value lambda*(sqrt(T+1))^dim is used, dim={1,2}.
%       accuracy [{'abs','rel'},'abs'] : absolute/relative accuracy mode.
%           We recommend accuracy='rel' with eps=0.1 (see below) if the
%           filter bandwidth is big enough (starting from L=T=6).
%       *eps [{>=0}] : precision up to which the problem must be solved,
%           absolute or relative depending on <accuracy>. We recommend 0.1
%           if accuracy='rel' and sigm*rho/2, where sigm is the noise
%           level, if accuracy='abs'.
%       tol [{>=0},1e-8] : numerical tolerance for comparisons. In most
%           cases 1e-8 is enough. 1e-12 and higher is not recommended, some
%           numerical instability is possible.
%       verbose [{0,1},0] : print output every <verbose> steps, never if 0.
%       solver [{'cvx','mp','nes','sub'},'nes'|'mp'] : solver. CVX, Mirror
%           Prox, Nesterov, or Projected Subgradient. 'nes' is assigned
%           by default if p=2, (it can be used only in that case), and 'mp'
%           otherwise. 'sub' can be used only if p=inf and not recommended.
%       squared [{0,1},0] : must be given only if solver='cvx',
%           constrained=0, and p=2. squared=0 corresponds to the penalized 
%           l_2, and squared=1 to the penalized l_2^2 minimization. Useful
%           to output performance curves for Nesterov (with squared=1) and
%           Mirror Prox (with squared=0), see TEST_MANUAL.
%       max_iter [{int},100000] : max number of iterations to perform.
%       max_cpu [{int},100000] : max computation time in seconds.
%       online [{0,1},1] : adaptive stepsize heuristic. Used if <solver> is
%           'mp' or 'nes'.
%       warm_start [{0,1},0] : initialize on signal DFT subspace solving a
%           linear regression problem by CVX. May be slow, not recommended.
%       phi0 : initial approximation of the filter DFT, normalized (i.e.
%           after change of variables). Zeros by default.
%       u0 : intial approximation of the dual variable when solver='mp',
%           normalized (i.e. after change of variables). Zeros by default.
%       l2_prox [{0,1},1] : use the Euclidean prox everywhere.
%       last_part [{int},50] : heuristic, average over the last iterations
%           (only if <solver> is 'mp').
%       if_trace_filter[{0,1},0] : return the matrix whose i-th column is
%           the filter computed at iteration i (only if <solver> is 'mp' or
%           'nes'). Useful to trace the true recovery error. Only for 1-D.
%       
%       
% Output params
%   sol : structure with the following fields
%       status : 1 if max_iter reeached, 2 if max_cpu reached, 0 otherwise.
%       cpu : total running time in sec
%       iter : total number of iterations done if <solver> is 'mp', 'nes', 
%           or 'sub'.
%       online_iter_rate : adaptive stepsize heuristic efficiency -- ratio
%           of internal and external interation counts.
%       online_stepsize_ratio :  adaptive stepsize heuristic efficiency --
%           ratio of theoretical and actual (average) stepsize. Only if
%           solver='mp'.
%       phi : reported primal solution, normalized (i.e. after the change
%           of variables).
%       u : reported dual solution, normalized (i.e. after the change of
%           variables). Only if solver='mp'.
%       pobj : primal objective value at the reported solution.
%       dobj : dual objective value at the reported solution. Only if
%           solver='mp'.
%       lowerb : lower bound on the minimum value. Only if solver='nes'.
%       gap : pobj-dobj (if solver='mp') or pobj-lowerb (if solver='nes').
%       filter : fitted filter, rho / sqrt(T+1)^(dim) * idft(sol.phi).
%       Filter : filters after each iteration, only if_trace_filter == 1.
%       Pobj : stats, primal objective value after each iteration, see
%           sol.pobj.
%       Dobj : stats, dual objective value after each iteration (if
%           solver='mp'), see sol.dobj.
%       Lowerb : stats, lower bound on the minimal value after each
%           iteration (if solver='nes'), see sol.lowerb.
%       Gap : stats, gap after each iteration (see sol.gap).
%       RelErrUb : stats, upper bound on the relative error after each
%           iteration.
%       Cpu : stats, time spent (in seconds) after each iteration.
%       Int_iter : stats, total number of internal iterations (due to the
%           adaptive stepsize search) after each iteration.
%
% Copyright : Dmitry Ostrovsky, 2017.
    
    % control.p, control.solver, and control.verbose are handled here
    if isfield(control,'p')
        if control.p ~= 2 && control.p ~= +inf
            error('p must be 2 or inf');
        end
    else
        control.p = 2; % because l_inf fitting is more conservative
    end
    if ~isfield(control,'solver')
        control.solver = 'mp';
        if control.p == 2
            control.solver = 'nes';
        end        
    end
    if ~isfield(control,'verbose')
        control.verbose = 0;
    end
    if strcmp(control.solver,'cvx')
        sol = fit_filter_cvx(y1,y2,control);
    elseif strcmp(control.solver,'sub')
        sol = fit_filter_subgrad(y1,y2,control);
    elseif strcmp(control.solver,'nes')
        sol = fit_filter_nesterov(y1,y2,control);
    elseif strcmp(control.solver,'mp')
        sol = fit_filter_mirror_prox(y1,y2,control);
    else
        error('Solver unknown');
    end
    if control.verbose
        fprintf('Summary\n');
        if ~sol.status
            fprintf('\tProblem status\t\t:   Solved\n');
        else
            fprintf('\tProblem status\t\t:   Unsolved - ');
            if sol.status == 1
                fprintf('max number of iterations reached\n');
            elseif sol.status == 2
                fprintf('max cpu time reached\n');
            end
        end
        fprintf('\tCPU Time\t\t:   %0.2f\n', sol.cpu);
        fprintf('\tIterations done\t\t:   %d\n', sol.iter);
        fprintf('\tPrimal objective\t:   %0.4f\n', sol.pobj);
        if strcmp(control.solver,'mp')
            fprintf('\tDual objective\t\t:   %0.4f\n', sol.dobj);
            fprintf('\tPrimal-dual gap\t\t:   %0.4f\n', sol.gap);
            fprintf('\tStepsize search rate\t:   %0.2f\n',... 
                sol.online_iter_rate);
            fprintf('\tStepsize ratio (th/pr)\t:   %0.2f\n',...
                sol.online_stepsize_ratio);
        end
    end
end