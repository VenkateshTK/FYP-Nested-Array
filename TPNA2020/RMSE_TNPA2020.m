% 2 level Nested array

% parameters
fc=2.4e9; % frequency of narrowband signal
c=physconst("lightspeed");
lambda=(c/fc); % wavelength
M1=9; % no.of elements in level1
M2=7; % no.of elements in level2
M=M2+M1; % total elemnts in subarray1
d=lambda/2; % distance between subarrays
d1=lambda/2; % separation in level 1
d2=(M1+1)*d1; % seperation in level 2
K=1; % no.of sources
sig_angle=[60,30]; % Angle of arrival (theta_1,phi_1;...;theta_k,phi_k)
% SNR=20; % SNR of signal
fs=20*fc; % sampling frequency
N=200; % no.of snapshots
t=0:1/fs:(N-1)*(1/fs);

% signal model construction
S=exp(1i*(2*pi*fc*t+(0:(pi/K):(K-1)*(pi/K))'));
% S=repmat(exp(1i*2*pi*fc*t),K,1);

A1=getManifoldMatrixA1(sig_angle(:,1),sig_angle(:,2),M1,M2,d1,lambda);
A2=A1*getPhi(sig_angle(:,1),sig_angle(:,2),d,lambda);


iterations=100;
SNR=0:20;
rmse=zeros(1,length(SNR));
for s=1:length(SNR)
    err=zeros(1,iterations);
    for i=1:iterations
        
        sigmansq=10^(-SNR(s)/10);
        n1=sqrt(sigmansq)*(randn(M,N)+1i*randn(M,N));
        n2=sqrt(sigmansq)*(randn(M,N)+1i*randn(M,N));
        
        x1=A1*S + n1;
        x2=A2*S + n2;
        
        % Proposed Algorithm
        
        % covariance and cross covariance matrices
        R11=(x1*x1')/N;
        R21=(x2*x1')/N;
        
        % Vectorizing (difference coarray)
        z1=R11(:);
        z2=R21(:);
        
        % forming virtual array by removing redundant elements
        idx=(getIndexOfUniqueElements(M1,M2));
        z1b=z1(idx);
        z2b=z2(idx);
        
        % toeplitz matrix
        Mb=M2*(M1+1);
        R11b=toeplitz(z1b(Mb:2*Mb-1),flip(z1b(1:Mb)));
        R21b=toeplitz(z2b(Mb:2*Mb-1),flip(z2b(1:Mb)));
        
        % Agumented covaraince matrix
        R=[R11b,R21b';R21b,R11b];
        
        % spectral function formulation
        [eigvec,~]=eigs(R,2*Mb);
        E=eigvec(:,K+1:2*Mb);
        
        % Spectral Estimation
        th=0:90;
        ph=270:360;
        P=zeros(length(th),length(ph));
        max=0;
        for m=1:length(th)
            for n=1:length(ph)
                ba=exp( 1i*(2*pi*d1/lambda)*(0:Mb-1)'*( sind(th(m))*sind(ph(n)) ) );
                bab=ba*exp(1i*(2*pi*d/lambda)*(cosd(th(m))*sind(ph(n))));
                b=[bab;ba];
                P(m,n)=10*log10(abs(1/(b'*(E*E')*b)));
                if (P(m,n)>max)
                    max=P(m,n);
                    th_k=th(m);
                    ph_k=360-ph(n);
                end
            end
        end
%         th_k
%         ph_k
        err(i)=mean(abs([th_k-sig_angle(:,1) ph_k-sig_angle(:,2)]));
%         surf(360-ph,th,(P),'EdgeColor',"none");
    end
    rmse(s)=rms(err);
end
figure;
semilogy(SNR,rmse)
xlabel('SNR in dB')
ylabel('RMSE')
title('RMSE for TPNA')
grid on

function uidx=getIndexOfUniqueElements(M1,M2)
narray=[0:M1 ((2:M2).*(M1+1)-1)];
cr=[];
for i=1:length(narray)
    cr=[cr -narray+narray(i)];
end
[~,uidx,~]=unique(cr);
end

function A=getManifoldMatrixA1(theta,phi,M1,M2,d1,lambda)
A=exp(1i*(2*pi*d1/lambda)*[0:M1 ((2:M2).*(M1+1)-1)]'*(sind(theta').*sind(phi')));
end

function Phi=getPhi(theta,phi,d,lambda)
Phi=diag(exp(1i*(2*pi*d/lambda)*(cosd(theta').*sind(phi'))));
end