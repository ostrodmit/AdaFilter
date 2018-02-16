function z_hat = l1_lr_prox_constrained_bisection(z, r, gamma, tol)
z = gamma * z;
y = abs(z);
if F(0,y,r) <= 0
    lambda = 0;
else % bisection
    lb = 0;
    rb = max(y);
    err = F(lb,y,r) - F(rb,y,r);
    iter = 0;
    max_iter = 100;
    while err > tol && iter <= max_iter
        iter = iter+1;
        lambda = (lb + rb)/2;
        if F(lambda,y,r) <= 0
            rb = lambda;
        else
            lb = lambda;
        end
        err = F(lb,y,r)-F(rb,y,r);
    end
end
z_hat = gen_soft_threshold(z, lambda, 1/(r-1));
end

function f = F(lambda,y,r)
    f = norm ( gen_soft_threshold(y,lambda,1/(r-1)), 1 ) - 1;
end