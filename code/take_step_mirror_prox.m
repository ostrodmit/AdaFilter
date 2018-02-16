function [phi, u, pphi, uu, nabla_pphi, nabla_uu,...
    prox_grad_phi_ext, prox_grad_phi_old] ... 
    = take_step_mirror_prox(constrained, phi, u, y1, b, rho, r, gamma,...
    pic, st, p, tol, lambda)
    r_phi = r(1); r_u = r(2); gamma_phi = gamma(1); gamma_u = gamma(2);
    
    vec_phi = phi(:); vec_u = u(:);
    prox_grad_phi = prox_grad(constrained, vec_phi, r_phi, gamma_phi);
    % save for the online stepsize check
    if ~pic
        prox_grad_phi_old = prox_grad_phi;
    else
        prox_grad_phi_old = vec2mat(prox_grad_phi);
    end
    prox_grad_u = prox_grad(1, vec_u, r_u, gamma_u);
    
    % evaluate gradients in the current point
    nabla_phi = adjoint_operator(u, rho, y1);
    nabla_u = direct_operator(phi, rho, y1) - b;

    % perform the intermediate step
    phiplus = prox_grad_phi - st * nabla_phi(:);
    pphi = l1_lr_prox(constrained, phiplus, r_phi,gamma_phi,tol,st,lambda);
%     % debug: check that the step was correct:
%     disp(norm(sign(pphi(:)) .* abs(pphi(:)).^(r_phi-1)...
%          -sign(vec_phi) .* abs(vec_phi).^(r_phi-1)));
%     disp(norm(gamma_phi * st * nabla_phi(:)))
    uplus = prox_grad_u + st * nabla_u(:);
    if p == 2 % project onto the unit l2-ball
        if norm(uplus,2) > 1
            uu = uplus / norm(uplus,2);
        else
            uu = uplus;
        end
    else
        uu = l1_lr_prox(1, uplus, r_u, gamma_u, tol);
    end
%     % debug: check that the step was correct:
%     disp(norm(sign(uu(:)) .* abs(uu(:)).^(r_u-1)...
%          -sign(vec_u) .* abs(vec_u).^(r_u-1)));
%     disp(norm(gamma_u * st * nabla_u(:)))
    if pic
        pphi = vec2mat(pphi); uu = vec2mat(uu);
    end
    
    % evaluate gradients at the intermediate point
    nabla_pphi = adjoint_operator(uu, rho, y1);
    nabla_uu = direct_operator(pphi, rho, y1) - b;
    
    prox_grad_pphi = prox_grad(constrained, pphi(:), r_phi, gamma_phi);
    % save for the online stepsize check
    if ~pic
        prox_grad_phi_ext = prox_grad_pphi;
    else
        prox_grad_phi_ext = vec2mat(prox_grad_pphi);
    end
    
    % perform the ultimate step
    if constrained
        phiplus = prox_grad_phi - st * nabla_pphi(:);
        phi = l1_lr_prox(1, phiplus, r_phi, gamma_phi, tol);
    else
        phiplus = prox_grad_pphi + st * nabla_phi(:) - st * nabla_pphi(:);
        % diff = st * (subgradient of the non-smooth component in phi_ext),
        % as follows from the optimality conditions at pphi
        phi = l1_lr_prox(0, phiplus, r_phi, gamma_phi, tol, st, 0);
    end
    uplus = prox_grad_u + st * nabla_uu(:);
    if p == 2
        if norm(uplus,2) > 1
            u = uplus/norm(uplus,2);
        else
            u = uplus;
        end
    else
        u = l1_lr_prox(1, uplus, r_u, gamma_u, tol);
    end
    if pic
        phi = vec2mat(phi); u = vec2mat(u);
    end
end