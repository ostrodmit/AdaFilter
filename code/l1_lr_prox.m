function z_hat = l1_lr_prox(constrained, z, r, gamma, tol, st, lambda)
if constrained
    z_hat = l1_lr_prox_constrained_newton_2(z, r, gamma, tol);
    % z_hat = l1_lr_prox_constrained_newton(z, r, gamma, tol); % jitter!
    % z_hat = l1_lr_prox_constrained_bisection(z, r, gamma, tol);
    % z_hat = l1_l2_prox_constrained_sort(z);% only if [r,gamma]==[2,1]
else
    z_hat = l1_lr_prox_penalized(z, r, gamma, st, lambda);
end