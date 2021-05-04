clc;
clear;

fc = 2.4e9;
c = physconst("lightspeed");
lambda = c/fc;
d = lambda/2;
M= 50;
obj = phased.ULA('NumElements',M,'ElementSpacing',d);
fs = 10*fc;
t=0:1/fs:4999*(1/fs);
S=[cos(2*pi*1e6*t').*cos(2*pi*fc*t') cos(2*pi*2e6*t').*cos(2*pi*fc*t') cos(2*pi*3e6*t').*sin(2*pi*fc*t)'];
D=size(S,2)
X=collectPlaneWave(obj,S,[10 30 70],fc);
sigmansq = 0.01;
N=sqrt(sigmansq)*(randn(length(t),M)+1i*randn(length(t),M));
X=X+N;
Rx=(X'*X)/length(t);
[eigvec,eigVal]=eigs(Rx,M);
eigv=eigvec(:,D+1:M);
syms theta
psi=(2*pi*d/lambda)*(sin(theta));
k=0:M-1;
a=exp(-1i*k'*psi);
p=1/((a'*eigv)*(eigv'*a));
th=(0:100);
theta=th*pi/180;
P=(double(subs(p)));
semilogy(th,abs(P));