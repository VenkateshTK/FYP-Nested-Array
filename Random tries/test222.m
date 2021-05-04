clc;
clear;

fc=2.4e9;
c=physconst("lightspeed");
lambda=(c/fc);
M1=3;
M2=3;
M=M2+M1;
d=lambda/2; %d = dx
d1=lambda/2; %d1 =dy
d2=(M1+1)*d1;
SNR=30; % in db
fs=10*fc;
t=0:1/fs:4999*(1/fs);

%signal  
signal_angle=[10 40 ;15 25 ];
freq = [pi/4 pi/3]';
S = exp(i*freq*t);
K=size(S,2);%no. of signals
sigmansq=10^(-SNR/10);
n1=sqrt(sigmansq)*(randn(length(t),M)+1i*randn(length(t),M));
n2=sqrt(sigmansq)*(randn(length(t),M)+1i*randn(length(t),M));



A1=getManifoldMatrixA1(signal_angle,M1,M2,d1,lambda);
Phi=getPhi(signal_angle,d,lambda);
A2=A1*Phi;

x1=S*(A1')+n1;
x2=S*(A2')+n2;

R11=cov(x1,1);
nv=eigs(R11,M);
sigmanest= mean(nv(K+1:M));
R11 = R11 - sigmanest*eye(M);

R21=(x2'*x1)/length(t);

Mb=M2*(M1+1);

z1 =[];
z2 =[];

d = [1,2,3,4,8,12];
m = [];
m.d1 = [];
m.d2 =[];
count = zeros(max(d),1);
for i = 1:length(d)
    for j = i: length(d)
        index = abs(d(i)-d(j));
        count(index+1) = count(index+1)+1;
        m(index+1).d1(count(index+1)) = i;
        m(index+1).d2(count(index+1)) = j;
    end
end

% -(Mb-1)
for i = 1:Mb-1 
    z1(i) =0;
    z2(i) =0;
    for j = 1:length(m(Mb-i+1).d2)
        z1(i) = z1(i) + R11(m(Mb-i+1).d2(j),m(Mb-i+1).d1(j));
        z2(i) = z2(i) + R21(m(Mb-i+1).d2(j),m(Mb-i+1).d1(j));
    end
    z1(i) = z1(i)/length(m(Mb-i+1).d2);
    z2(i) = z2(i)/length(m(Mb-i+1).d2);
end 

for i = Mb:1:2*Mb-1 
    z1(i) =0;
    z2(i) =0;
    for j = 1:length(m(i-Mb+1).d1)
        z1(i) = z1(i) + R11(m(i-Mb+1).d1(j),m(i-Mb+1).d2(j));
        z2(i) = z2(i) + R21(m(i-Mb+1).d1(j),m(i-Mb+1).d2(j));
    end
    z1(i) = z1(i)/length(m(i-Mb+1).d2);
    z2(i) = z2(i)/length(m(i-Mb+1).d2);
end 

z1b = z1';
z2b = z2';

Mb=M2*(M1+1);

c1 = z1b(1:Mb);
r1 = z1b(Mb:length(z1b));

Va = (hankel(c1,r1));

c2 = z2b(1:Mb);
r2 = z2b(Mb:length(z2b));

Vb = (hankel(c2,r2));

D = Vb*pinv(Va);

[F,psi] = eigs(D,Mb);

u = getUvector(F,Mb,K,d1,lambda);
v = getVvector(psi,K,d1,lambda);

ThetaMat = asin(sqrt(u.^2 + v.^2))*180/pi
PhiMat = angle(u+v*1i)*180/pi


function U = getUvector(F,Mb,K,d,lambda)
    dr = 2*pi*d/lambda;
    U = [];
    for k=1:K
       U(k) = (1/((Mb-1)*dr))*(sum(angle(F(2:Mb,k)./F(1:Mb-1,k))));
    end
end

function V = getVvector(psi,K,d1,lambda)
    dr = 2*pi*d1/lambda;
    V = [];
    for k=1:K
        V(k) = angle(psi(k,k))/dr;
    end
end


function A=getManifoldMatrixA1(angle,M1,M2,d1,lambda)
A=exp(1i*(2*pi*d1/lambda)*[0:M1 ((2:M2).*(M1+1)-1)]'*(sind(angle(1,:)).*sind(angle(2,:))));
end

function Phi=getPhi(angle,d,lambda)
Phi=diag(exp(1i*(2*pi*d/lambda)*(cosd(angle(1,:)).*sind(angle(2,:)))));
end