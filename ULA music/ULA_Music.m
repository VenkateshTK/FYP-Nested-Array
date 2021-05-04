fc=2.4e9;
c=299792458;
lambda=(c/fc);
d=lambda/2;
M=30;
SNR=20;
A=getULASteeringVec([10 30 60].*pi/180,M,d,lambda);
fs=10*fc;
t=0:1/fs:4999*(1/fs);
S=[cos(2*pi*1e6*t').*cos(2*pi*fc*t') cos(2*pi*2e6*t').*cos(2*pi*fc*t') cos(2*pi*3e6*t').*cos(2*pi*fc*t)'];
D=size(S,2);
sigmansq=10^(-SNR/10);
N=sqrt(sigmansq)*(randn(length(t),M)+1i*randn(length(t),M));
X=S*A'+N;
Rx=(X'*X)/length(t);
[eigvec,eigVal]=eigs(Rx,M);
eigv=eigvec(:,D+1:M);
th=(0:90);
P=zeros(1,length(th));
for i=1:length(th)
    a=getULASteeringVec(th(i)*pi/180,M,d,lambda);
    P(i)=abs(1/((a'*eigv)*(eigv'*a)));
end
plot(th,10*log10(P));
title('Angle of Arrival estimate using MUSIC')
xlabel('\theta in degrees')
ylabel('|P(\theta)| in dB')
grid on

function A=getULASteeringVec(theta,num_elements,dist,lambda)
A=exp(-1i*(2*pi*dist/lambda)*(0:num_elements-1)'*sin(theta));
end
