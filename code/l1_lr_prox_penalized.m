function z_hat = l1_lr_prox_penalized(z, r, gamma, st, lambda)
n = length(z);
s = r/(r-1);
C = n^((2-r)/s)/gamma;
y = subplus(abs(z)-st*lambda);
if r == 2 && gamma == 1,
    z_hat = sign(z) .* y;
else
    norm_y_s = norm(y,s);
    if  norm_y_s > 0
        z_hat = sign(z) .* (y * norm_y_s^(r-2)).^(1/(r-1)) / C;
    else
        z_hat = zeros(size(z));
    end
end
end