% 2 level Nested array
fc=2.4e9;
c=physconst("lightspeed");
lambda=(c/fc);
M1=2;
M2=2;
M=M2+M1;
d=lambda/2; %d = dx
d1=lambda/2; %d1 =dy
d2=(M1+1)*d1;
SNR=20;
fs=10*fc;
t=0:1/fs:4999*(1/fs);
S=[cos(2*pi*2e6*t').*cos(2*pi*fc*t') cos(2*pi*3e6*t').*cos(2*pi*fc*t)'];
K=size(S,2);%no. of signals
sigmansq=10^(-SNR/10);
n1=sqrt(sigmansq)*(randn(length(t),M)+1i*randn(length(t),M));
n2=sqrt(sigmansq)*(randn(length(t),M)+1i*randn(length(t),M));
signal_angle=[40 80; 10 30];
A1=getManifoldMatrixA1(signal_angle,M1,M2,d1,lambda);
Phi=getPhi(signal_angle,d,lambda);
A2=A1*Phi;
x1=S*A1'+n1;
x2=S*A2'+n2;

R11=cov(x1,1);
[~,nv]=eigs(R11,M);
sigmanest=mean(nv(K+1:M));
R11 = R11 - sigmanest*eye(M);
z1=R11(:);
z1b=flip(z1(getIndexOfUniqueElements(M1,M2))); %zib = va bar
R21=(x2'*x1)/length(t);


z2=R21(:);
z2b=flip(z2(getIndexOfUniqueElements(M1,M2))); %z2b = vb bar
Mb=M2*(M1+1);

c = z1b(1:Mb);
r = z1b(Mb:length(z1b));

Va = hankel(c,r);

c = z2b(1:Mb);
r = z2b(Mb:length(z2b));

Vb = hankel(c,r);

D = Vb*pinv(Va);

[F,psi] = eigs(D,Mb);

u = getUvector(F,Mb,K,d,lambda)
v = getVvector(psi,K,d1,lambda)

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

function uidx=getIndexOfUniqueElements(M1,M2)
narray=[0:M1 ((2:M2).*(M1+1)-1)];
cr=[];
for i=1:length(narray)
    cr=[cr -narray+narray(i)];
end
[~,uidx,~]=unique(cr);
end
function A=getManifoldMatrixA1(angle,M1,M2,d1,lambda)
A=exp(1i*(2*pi*d1/lambda)*[0:M1 ((2:M2).*(M1+1)-1)]'*(sind(angle(1,:)).*sind(angle(2,:))));
end
function Phi=getPhi(angle,d,lambda)
Phi=diag(exp(1i*(2*pi*d/lambda)*(cosd(angle(1,:)).*sind(angle(2,:)))));
end