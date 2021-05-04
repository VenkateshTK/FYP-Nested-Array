% 2 level Nested array
clc;
clear;

fc=2.4e9;
c=physconst("lightspeed");
lambda=(c/fc);
M1=3;
M2=3;
M=M2+M1;
d=lambda/2;
d1=lambda/2;
d2=(M1+1)*d1;
SNR=20;

fs=10*fc;
t=0:1/fs:(1e5-1)*(1/fs);
S=[cos(2*pi*1e6*t').*cos(2*pi*fc*t') cos(2*pi*2e6*t').*cos(2*pi*fc*t')];
D=size(S,2);

sigmansq=10^(-SNR/10);
n1=sqrt(sigmansq)*(randn(length(t),M)+1i*randn(length(t),M));
n2=sqrt(sigmansq)*(randn(length(t),M)+1i*randn(length(t),M));

angle=[70 60; 50 40]
A1=getManifoldMatrixA1(angle,M1,M2,d1,lambda);
Phi=getPhi(angle,d,lambda);
A2=A1*Phi;

x1=S*A1'+n1;
x2=S*A2'+n2;

R11= x1'*x1/length(t);
z1=R11(:);
index_array = flip(getIndexOfUniqueElements(M1,M2))
z1b=z1(index_array);

R21=(x2'*x1)/length(t);
z2=R21(:);
z2b=z2(index_array);

Mb=M2*(M1+1);
c=z1b(Mb:2*Mb-1);
r=flip(z1b(1:Mb));
R11b=toeplitz(c,r);

c=z2b(Mb:2*Mb-1);
r=flip(z2b(1:Mb));
R21b=toeplitz(c,r);
R=[R11b R21b'; R21b R11b];

[eigvec,eigVal]=eigs(R,2*Mb);
E=eigvec(:,D+1:2*Mb);

th=0:90;
ph=0:90;
P=[];

for m=1:length(th)
    for n=1:length(ph)
        ba=exp(1i*(2*pi*d1/lambda)*(0:Mb-1)'*(sind(th(m))*sind(ph(n))));
        bab=ba*exp(1i*(2*pi*d/lambda)*(cosd(th(m))*sind(ph(n))));
        b=[ba;bab];
        P(m,n)=10*log10((abs(1/(b'*(E*E')*b))));
    end
end

figure;
surf(th,ph,P,'EdgeColor',"none");
title("Power Spectrum in TPNA, M1=3, M2=3");
xlabel("theta");
ylabel("phi");
zlabel("Spectrum");
% alph=acosd(sin(angle(1,:)).*sin(angle(2,:)))
% be=acosd(cos(angle(1,:)).*sin(angle(2,:)))

En1=E(1:Mb,:);
En2=E(Mb+1:2*Mb,:);
 
v=sym('v');
b1=power(v,[0:Mb-1]');
p1=(b1')*(En1*(En1'))*b1;
p2=(b1')*(En1*(En2'))*b1;
p3=(b1')*(En2*(En1'))*b1;
p4=(b1')*(En2*(En2'))*b1;
det_q=p1*p4-p2*p3;
%[vest,param,cond]=solve(det_q==0,v,"ReturnConditions",true)


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