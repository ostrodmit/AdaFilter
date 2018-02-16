function scale=scale_block_l1(rdata,Eucl)
scale.n=rdata.n;
if Eucl,
    scale.p=2;
    scale.cf=1;
else
    n=rdata.n;
    if n<=2,
        if n==1,
            gamma=1;
        else
            gamma=0.5;
        end;
        scale.p=2;
    else
        scale.p=1+1/log(n);
        gamma=1/(exp(1)*log(n));
    end;
    deg=(scale.p-1)*(2-scale.p)/scale.p;
    scale.cf=(n^deg)/gamma;
end;
