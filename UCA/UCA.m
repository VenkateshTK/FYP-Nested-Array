clc;
clear;

fc = 2.4e9;
c = physconst("lightspeed");
lambda = c/fc;
M = 20;
r = M*lambda/(4*pi);
sUCA = phased.UCA(M,r);
viewArray(sUCA);

fs = 10*fc;
t=0:1/fs:3999*(1/fs);
S=[cos(2*pi*1e6*t').*cos(2*pi*fc*t') cos(2*pi*2e6*t').*cos(2*pi*fc*t') cos(2*pi*3e6*t').*sin(2*pi*fc*t)'];
D=size(S,2)

X = collectPlaneWave(sUCA,S,[[10;10] [40;30] [80;40]],fc);
sigmansq = 0.01;
N=sqrt(sigmansq)*(randn(length(t),M)+1i*randn(length(t),M));
X=X+N;

Rx=(X'*X)/length(t);
[eigvec,eigVal]=eigs(Rx,M);
eigv=eigvec(:,D+1:M);

p=zeros(90);
th =(1:0.5:90);
ph =(1:0.5:90);
for i = 1:length(th)
    for j = 1:length(ph)
        k=1:2:2*M;
        theta=th(i)*pi/180;
        phi = ph(j)*pi/180;
        phase = (2*pi*r/lambda)*(sin(theta)*cos(phi-(pi*(k')/M)));
        a=exp(1i*phase);
        p(i,j)=10*log(abs(1/((a'*eigv)*(eigv'*a))));
    end
end
figure;
surf(th,ph,p,'EdgeColor',"none");


estimator = phased.MUSICEstimator2D('SensorArray',sUCA,...
    'OperatingFrequency',fc,...
    'NumSignalsSource','Property',...
    'DOAOutputPort',true,'NumSignals',3,...
    'AzimuthScanAngles',0:.5:90,...
    'ElevationScanAngles',0:.5:90);
[~,doas] = estimator(X)

figure;
plotSpectrum(estimator);
