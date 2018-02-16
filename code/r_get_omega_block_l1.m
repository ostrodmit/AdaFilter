function [omega,domega]=r_get_omega_block_l1(arg,scale)
%
aarg=abs(arg);
mx=max(aarg);
if mx==0,
    omega=0;
    domega=0*arg;
    return;
end;
tmp=sum((aarg/mx).^scale.p);
domega=zeros(scale.n,1);
fctr=mx*(tmp^(2/scale.p-1));
nind= aarg>0;
%domega(nind)=fctr*arg(nind)./aarg(nind).*(aarg(nind)/mx).^(scale.p-1);
domega(nind)=fctr*sign(arg(nind)).*(aarg(nind)/mx).^(scale.p-1);
omega=scale.cf*arg'*domega/2;
domega=scale.cf*domega;
end % endof GetOmegaBlockEll1