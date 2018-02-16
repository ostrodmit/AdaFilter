function [y,f,ft]=ysim(t,Snr,type)
% generate signal of length t
% noise variance is defined by the SNR 
xx=[1:t]';
tt=[1:10*t]';
no=nargout;
if no>2
    xt=tt/10;
end

if type==1,
    f=sin([1:t]*sqrt(200))+sqrt(5)*sin([1:t]*sqrt(337));
    f=f(:);
elseif type==2,
    f=sin([1:t]/sqrt(17))';
elseif type==3,
    f=4*sqrt(t-xx)/sqrt(t);
elseif type==4,
    f=2*(2/t*xx-1.8*xx.^2/t^2).*sin(xx'*sqrt(55))';
elseif type==5,
    f=4*(2/t*xx-3.2*xx.^2/t^2+1.5*xx.^4/t^4).*sin(xx'*sqrt(45))';
elseif type==6,
    f=4*(2/t*xx-3*xx.^2/t^2+1.2*xx.^4/t^4).*sin(9*xx'/t)';
elseif type==7,
    f=-0.7191*cos(0.0017*xx)+1.5200*cos(0.0025*xx)...
        +0.6041*cos(0.3526*xx)-0.0676*cos(0.8303*xx)-0.1892*cos(3.2371*xx);
    if no>2,
        ft=-0.7191*cos(0.0017*xt)+1.5200*cos(0.0025*xt)...
            +0.6041*cos(0.3526*xt)-0.0676*cos(0.8303*xt)-0.1892*cos(3.2371*xt);
    end
elseif type==8,
    gt=1+cumsum(1/2/sqrt(10*t)*randn(t*10,1));
    g=gt(10:10:10*t);
    %    f=-4*(2/t*xx-3.2*xx.^2/t^2+1.5*xx.^4/t^4).*(0.7191*cos(0.0017*xx)+...
    %        1.5200*cos(0.0025*xx)+0.6041*cos(0.3526*xx)-0.0676*cos(0.8303*xx)-0.1892*cos(3.2371*xx));
    f=g(xx).*(-0.7191*cos(0.0017*xx)+...
        1.5200*cos(0.0025*xx)+0.6041*cos(0.3526*xx)-0.0676*cos(0.8303*xx)-0.1892*cos(3.2371*xx));
    if no>2,
        ft=gt.*(-0.7191*cos(0.0017*xt)+...
            1.5200*cos(0.0025*xt)+0.6041*cos(0.3526*xt)-0.0676*cos(0.8303*xt)-0.1892*cos(3.2371*xt));
    end
end
sig=norm(f)/sqrt(t)/Snr;
y=f+sig*randn(t,1);


%