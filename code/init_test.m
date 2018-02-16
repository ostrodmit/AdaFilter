T = 50; L = 49;
k = 10; rho = 20;
sigma = 1;
%generating a 2d signal comprised of k sines
freqsX = rand(k,1);
freqsY = rand(k,1);
n = T+L+1;
x=0;
for i = 1:k
    a = repmat(2*pi*1i*freqsX(i)*(1:n), n, 1);
    b = repmat(2*pi*1i*freqsY(i)*(1:n)', 1, n);
    x = x+randn(1,1)*exp(a+b);
end
y = x + sigma * randn(n,n);

[n,m]=size(y);
imsize=n*m;
tic
z2=fft2(y)/sqrt(imsize);

z1=reshape(z2,imsize,1);
[sz,iz]=sort(abs(z1), 'descend');
ssz=cumsum(sz); 
issz=floor(rho);
ez1=zeros(imsize,1); 
ez1(iz(1:issz))=z1(iz(1:issz));
toc
disp(['initial inf-norm: ', num2str(max(abs(z1))),'; new inf-norm: ', num2str(max(abs(z1-ez1)))])

ez2=reshape(ez1,n,m);
subplot(2,1,1); imagesc(50*abs(z2))
title('Original FFT')
subplot(2,1,2); imagesc(50*abs(ez2))
title('Truncated FFT')

