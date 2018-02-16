format long
SNR=1;
L=1000;T=1000;
K=100;
%om=randn(2,1);
% s=sin(om(2)*[1:L+T-1]/sqrt(5));;s=s(:);
% y=s+randn(length(s),1)*norm(s)/SNR/sqrt(length(s)); y=y(:);

stype=7;
[~,s]=ysim(L+T+K+1,SNR,stype);
%Dima's change
sigm = norm(s)/sqrt(L+T+K+1)/SNR;
y = s + sigm * ones(L+T+K+1,1);
rho=10;
p = 2;
control = struct('p', p, 'rho', rho, 'solver', 'nes', 'constrained', 1,... 
    'accuracy', 'rel', 'eps', 0.1, ...
    'tol', 1e-8, 'max_iter', 1000000, 'max_cpu', 1000, 'l2_prox', 1,... 
    'warm_start', 0, 'online', 1, 'last_part', 0, 'verbose', 0);
tic
sol=fit_filter(y(1:L+T+1),y(T+1:L+T+1),control);
sol.iter
% control.ItrMax=1000;
% control.print=1000;
% control.gap=1/3;
toc
tic
for i=1:K
    control.phi0 = sol.phi;
    if strcmp(control.solver, 'mp') 
        control.u0 = sol.u;
    end
    sol=fit_filter(y(1+i:L+T+1+i),y(T+1+i:L+T+1+i),control);
    sol.iter
end
toc