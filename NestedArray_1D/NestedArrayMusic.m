fc=2.4e9; % frequency of narrowband signal
c=physconst("lightspeed");
lambda=(c/fc); % wavelength
M1=3; % no.of elements in level1
M2=2; % no.of elements in level2
M=M2+M1; % total elemnts in subarray1
d1=lambda/2; % separation in level 1
d2=(M1+1)*d1; % seperation in level 2
K=6; % no.of sources

sig_angle=[5;30;50;70;90;10]; % Angle of arrival (theta_1;...;theta_k)

SNR=20; % SNR of signal
fs=50*fc; % sampling frequency
N=5000; % no.of snapshots
t=0:1/fs:(N-1)*(1/fs);

% signal model construction
S=exp(1i*(2*pi*fc*t+(0:(pi/K):(K-1)*(pi/K))'));

A1=getManifoldMatrixA(sig_angle,M1,M2,d1,lambda);

sigmansq=10^(-SNR/10);
n1=sqrt(sigmansq)*(randn(M,N)+1i*randn(M,N));
x1=A1*S + n1;

% covariance and cross covariance matrices
R11=(x1*x1')/N;

z1=R11(:);
z1b=z1(getIndexOfUniqueElements(M1,M2));

Mb=M2*(M1+1);
R=toeplitz(z1b(Mb:2*Mb-1),flip(z1b(1:Mb)));

[eigvec,~]=eigs(R,Mb);
E=eigvec(:,K+1:Mb);

th=(0:180);
P=zeros(1,length(th));
for i=1:length(th)
    a=getVirtualULAsteeringVec(th(i),Mb,d1,lambda);
    P(i)=abs(1/((a'*E)*(E'*a)));
end
plot(th,10*log10(P));
title('MUSIC for nested array with M_1 = 3 M_2 = 2')
xlabel('\theta in degrees')
ylabel('|P(\theta)| in dB')
grid on

function uidx=getIndexOfUniqueElements(M1,M2)
narray=[0:M1 ((2:M2).*(M1+1)-1)];
cr=[];
for i=1:length(narray)
    cr=[cr -narray+narray(i)];
end
[~,uidx,~]=unique(cr);
end

function A=getVirtualULAsteeringVec(theta,Mb,d1,lambda)
A=exp(1i*(2*pi*d1/lambda)*(0:Mb-1)'*cosd(theta));
end

function A=getManifoldMatrixA(theta,M1,M2,d1,lambda)
A=exp(1i*(2*pi*d1/lambda)*[0:M1 ((2:M2).*(M1+1)-1)]'*(cosd(theta')));
end