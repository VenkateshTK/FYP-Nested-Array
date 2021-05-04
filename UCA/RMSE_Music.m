fc=2.4e9;
c=299792458;
lambda=(c/fc);
M=20;
r = M*lambda/(4*pi);

fs = 10*fc;
t=0:1/fs:4999*(1/fs);
S=[cos(2*pi*1e6*t').*cos(2*pi*fc*t') cos(2*pi*2e6*t').*cos(2*pi*fc*t')];
D=size(S,2);
SNR=0:2:20;
angle=[60 40;20 40];
rmse=zeros(1,length(SNR));
for j=1:length(SNR)
    trials=100;
    er=zeros(1,trials);
    for i=1:trials
        A=getUCASteeringVec(angle.*pi/180,M,r,lambda);
        sigmansq=10^(-SNR(j)/10);
        N=sqrt(sigmansq)*(randn(length(t),M)+1i*randn(length(t),M));
        X=S*A'+N;
        doa=musicUcaEstimator(X,M,r,lambda,D);
        error=sortrows(angle')-sortrows(doa');
        er(i)=norm(error(:))/power(length(error(:)),2);
    end
    rmse(j)=norm(er)/(trials^2);
end
semilogy(SNR,rmse)
xlabel('SNR in dB')
ylabel('RMSE')
title('RMSE in UCA')

function doa=musicUcaEstimator(X,M,r,lambda,D)
Rx=(X'*X)/size(X,1);
[eigvec,~]=eigs(Rx,M);
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

% surf(th,ph,p,'EdgeColor',"none");
[~,idx]=maxk(p(:),D);
doa(1,:)=th(mod(idx,length(ph))+1);
doa(2,:)=ph(ceil(idx/length(th)));
end
function A=getUCASteeringVec(angle,num_elements,r,lambda)
theta=angle(1,:);
phi=angle(2,:);
phase = (2*pi*r/lambda)*( sin(theta).*cos(phi-(pi*((0:num_elements-1)')/num_elements)) );
A=exp(1i*phase);
end