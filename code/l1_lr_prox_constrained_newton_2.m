function z_hat = l1_lr_prox_constrained_newton_2(z, p, gamma, tol)
% % \|x\|_p^2/2 prox
ni=nargin;
if ni<4, tol=1e-8; end
s=1/(p-1); q=p/(p-1);
% renormalize
z = gamma/length(z)^((2-p)/q) * z;
y = abs(z);
temp=max(y);
if temp==0, z_hat=zeros(size(z)); return; end
lambda = max(temp-1-1e-20,0);
ytemp=subplus(y-lambda);
if p==2
    ytemp2=ytemp; 
else
    ymax=max(ytemp); ytemp=ytemp/ymax;
    temp=sum(ytemp.^q);
    ytemp2=ymax*ytemp.^s/temp^(2/p-1); 
end
f = sum(ytemp2) - 1;
max_iter = 100; iter=0;
while f > tol && iter <= max_iter
    if p==2, fprim = -nnz(ytemp);
    else
        fprim = -((q-1)*sum(ytemp.^(q-2))*temp-(q-2)*sum(ytemp.^s)^2)/temp^(2/p);
    end
    lambda = lambda - f/fprim;
    ytemp = subplus(y-lambda);
    if p==2, ytemp2=ytemp; 
    else
        % stabilized computation
        ymax=max(ytemp); ytemp=ytemp/ymax;
        temp=sum(ytemp.^q); 
        ytemp2=ymax*ytemp.^s/temp^(2/p-1);
    end
    f = sum(ytemp2)-1;
    iter = iter+1;
end
z_hat = sign(z).*ytemp2;
%    iter
end %endof l1_lr_prox_constrained_newton_ai