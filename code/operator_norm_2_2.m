function opnorm = operator_norm_2_2(pic,y1,T,rho,nrep)
% computing the leading eigenvector of A^*A
phi = randn(T+1,(T+1)^pic) + 1i * randn(T+1,(T+1)^pic);
if nargin == 4, nrep = 100; end
for i = 1:nrep
    u = direct_operator(phi, rho, y1);
    phi = adjoint_operator(u, rho, y1);
    phi = phi / norm(phi(:),2);
end
u = direct_operator(phi, rho, y1);
phi = adjoint_operator(u, rho, y1);        
opnorm = sqrt(norm(phi(:),2)); % operator norm of A