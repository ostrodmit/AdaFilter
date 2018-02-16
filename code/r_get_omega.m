function [omega,domega]=r_get_omega(arg,scale)
mx=max(abs(arg));
if mx==0,
    omega=0;
    domega=0*arg;
    return;
end;
tmp=sum((abs(arg)/mx).^scale.p);
domega=mx*(tmp^(2/scale.p-1))*((abs(arg)/mx).^(scale.p-1)).*sign(arg);
omega=scale.cf*arg'*domega/2;
domega=-scale.cf*domega;
end % endof r_get_omega