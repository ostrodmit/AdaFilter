function nabla_omega = prox_grad(constrained, z, r, gamma)
% the gradient of the prox function
if r == 2 && gamma == 1
    nabla_omega = z;
else
    constrained = 0;
    if constrained
        nabla_omega = sign(z) .* abs(z).^(r-1) / gamma;
    else
        n = length(z);
        nabla_omega = norm(z,r)^(2-r) * sign(z).*abs(z).^(r-1) ...
            * n^((r-1)*(2-r)/r) / gamma;
    end
end
end