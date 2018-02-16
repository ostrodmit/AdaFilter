function [] = test_auto(scenario,snr,rho,algos)
test_input.sce = scenario;
if nargin > 1
    test_input.snr = snr;
end
if nargin > 2 
    test_input.rho = rho;
end
if nargin > 3
    test_input.alg = algos;
end 
test_manual(test_input);