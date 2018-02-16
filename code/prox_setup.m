function [r, gamma] = prox_setup(l2,n)
r = 2;
gamma = 1;
if ~l2
    if n == 2
        gamma = 0.5;
    end
    if n > 2
        r = 1 + 1/log(n);
        gamma = 1/(exp(1)*log(n));
    end
end