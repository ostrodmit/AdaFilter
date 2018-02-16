function zhat= l1_l2_prox_constrained_sort(z)
if sum(abs(z))<=1, zhat=z;
else
    % compute lambda
    n = length(z);
    y = sort(abs(z), 'descend');
    temp=cumsum(y)-(1:n)'.*y;
    tind=(temp<=1);m=sum(tind);
    lambda=(sum(y(tind))-1)/m;
    % now compute zhat
    zhat=sign(z).*subplus(abs(z)-lambda);
    
end
end%endof l1_l2_prox_constrained_sort

