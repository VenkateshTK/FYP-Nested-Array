fc=2.4e9;
c=299792458;
lambda=(c/fc);
M=30;
r = M*lambda/(4*pi);
SNR=20;
fs = 10*fc;
t=0:1/fs:4999*(1/fs);
S=[cos(2*pi*1e6*t').*cos(2*pi*fc*t') cos(2*pi*2e6*t').*cos(2*pi*fc*t') cos(2*pi*3e6*t').*sin(2*pi*fc*t)'];
D=size(S,2);
A=getUCASteeringVec([10 40 80; 10 30 40].*pi/180,M,r,lambda);
sigmansq=10^(-SNR/10);
N=sqrt(sigmansq)*(randn(length(t),M)+1i*randn(length(t),M));
X=S*A'+N;
Rx=(X'*X)/length(t);
[eigvec,eigVal]=eigs(Rx,M);
eigv=eigvec(:,D+1:M);

th =(1:0.5:90);
ph =(1:0.5:90);
p=zeros(length(th),length(ph));
for i = 1:length(th)
    for j = 1:length(ph)
        a=getUCASteeringVec([th(i); ph(j)].*pi/180,M,r,lambda);
        p(i,j)=10*log10(abs(1/((a'*eigv)*(eigv'*a))));
    end
end

figure;
surf(th,ph,p,'EdgeColor',"none");
title("Power Spectrum in UCA");
xlabel("theta");
ylabel("phi");
zlabel("Spectrum in db");

function A=getUCASteeringVec(angle,num_elements,r,lambda)
theta=angle(1,:);
phi=angle(2,:);
phase = (2*pi*r/lambda)*( sin(theta).*cos(phi-(pi*((0:num_elements-1)')/num_elements)) );
A=exp(1i*phase);
end
