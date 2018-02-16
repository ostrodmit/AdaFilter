function [z,omegaz,domegaz]=prox_block_l1(eta,alpha,lambda,scale)
scaleT=scale;
scaleT.cf=1;
scaleT.p=scale.p/(scale.p-1);
aarg=abs(eta);
arg=max(aarg-alpha*lambda,0)/scale.cf;

[~,w]=r_get_omega(arg,scaleT);
nind=arg>0;
z=zeros(scale.n,1);
z(nind)=w(nind)./aarg(nind).*eta(nind>0);
[omegaz,domegaz]=r_get_omega_block_l1(z,scale);
end % endof prox_blok_l1
