clc;
clear;

iterations = 100;
real_value = 10;
rmse = [];
it = 1;

fc = 2.4e9;
c = physconst("lightspeed");
lambda = c/fc;
d = lambda/2;
M = 6;
array = phased.ULA('NumElements',M,'ElementSpacing',d);
fs = 10*fc;
t=0:1/fs:4999*(1/fs);
S=[cos(2*pi*1e6*t').*cos(2*pi*fc*t')];
D=size(S,2);
X=collectPlaneWave(array,S,[real_value],fc);


for SNR = -20 : 1 : 2
    sum = 0;
    snr_nw = 10^(SNR/10);
    sigmansq = 1/sqrt(snr_nw);

    for i = 1: iterations
        
        N=sqrt(sigmansq)*(randn(length(t),M)+1i*randn(length(t),M));
        estimator = phased.MUSICEstimator('SensorArray',array,...
        'OperatingFrequency',fc,...
        'DOAOutputPort',true,'NumSignalsSource','Property',...
        'NumSignals',1);
        [y,doas] = estimator(X+N);
        est = doas(1);
        sum = sum + (real_value-est)^2;

    end

    rmse(it) = sqrt(sum/iterations);
    out = rmse(it)
    it = it+1;
end

figure;
semilogy((-20:1:2),rmse);
title(" RMSE in ULA");
xlabel("SNR in db");
ylabel("RMSE")
