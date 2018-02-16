X=randn(500,500);
G=gpuArray(X);
tic; for i=1:100, G1=fft2(G); end; toc
tic; for i=1:100, G1=fft2(G); end; X1=gather(G1); toc
tic; for i=1:100, X1=fft2(X); end; toc