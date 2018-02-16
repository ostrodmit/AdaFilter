function z_hat = l1_lr_prox_constrained_newton(z, r, gamma, tol)
z = gamma * z;
y = abs(z);
lambda = max(max(y)-1-1e-20,0);
q=1/(r-1); s=2-r;
ytemp=subplus(y-lambda);
if r~=2, ytemp=ytemp.^q; end
f = sum(ytemp) - 1;
max_iter = 100; iter=0;
while f > tol && iter <= max_iter
    iter = iter+1;
    if r==2, fprim = -nnz(ytemp); 
    else fprim = -sum(ytemp.^s)/(r-1); 
    end
    lambda = lambda - f/fprim;
    ytemp = subplus(y-lambda);
    if r~=2, ytemp = ytemp.^q; end
    f = sum(ytemp)-1;
end
%iter
z_hat = sign(z).*ytemp;