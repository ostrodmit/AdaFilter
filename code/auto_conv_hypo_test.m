clear all
close all
seed=2;
s = 1;
n = 10;
ntrials = 100;
rng(seed)
omega = 2*pi*rand(s,1);
%omega = 2*pi*(1:s)'/(n+1);
t = (0:n);
T = (0:2*n);
x = sum(exp(1i*omega*t),1);
xl = sum(exp(1i*omega*T),1);
%xx = conv(x,x);
xx = cconv(x,x,n+1);
figure
plot(real(fft(xl)))
figure
plot(real(fft(xx)))