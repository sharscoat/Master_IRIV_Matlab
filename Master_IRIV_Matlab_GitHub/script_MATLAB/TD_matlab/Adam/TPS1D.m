
clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Partie 1:

te=0.00001;
naq=10001;
taq=(naq-1)*te;
t =[0:te:taq];

Amp=1000
nu=1000
phi=0

Y=Amp*cos(2*pi*nu*t+phi);
figure
plot(Y)
plot(t,Y)


sY=fft(Y);
sY=fftshift(fft(Y));

RY=real(sY);
IY=imag(sY);
MY=abs(sY);

nue=1/te;
numax=nue/2;
dnu=1/taq;
nus=[-numax:dnu:numax]
plot(nus,MY)

% Y1=4*cos(2*pi*1000*t);
% Y2=2*cos(2*pi*2500*t+(pi/2));
% Y3=3*cos(2*pi*9000*t+pi);
% 
% Yt=Y1+Y2+Y3;
% 



Amp=[4 2 3 6];
nu = [1000 2500 9000 45000];
phi = [0 pi/2 pi 0]; 

Yt= zeros(1,naq);


for i = length(Amp)
    Yt= Yt+Amp(i)*cos(2*pi*nu(i)*t+phi(i))
end

sYt=fftshift(fft(Yt));
MYt=abs(sYt);
FYt=unwrap(angle(sYt));
RYt=real(sYt);
IYt=imag(sYt);
plot(nus,MYt);


figure;
subplot(2,2,1)
plot(t,Yt)
subplot(2,2,2)
plot(nus,MYt)
subplot(2,2,3)
plot(nus,RYt)
subplot(2,2,4)
plot(nus,FYt)


% for i = 1:(naq-1)/20
%     porte(i)=1
% end
% 
% for i=((naq-1)/20+1):naq
%     porte(i)=0
% end

Te=20
porte = zeros(1,naq);
porte(1:(naq-1)/Te)=1;


Y_2=porte.*Yt



sY_2=fftshift(fft(Y_2));
RY_2=real(sY_2);
IY_2=imag(sY_2);
MY_2=abs(sY_2);


 figure;
 subplot(2,2,1)
 plot(t,Y)
 subplot(2,2,2)
 plot(nus,MY_2)
 subplot(2,2,3)
 plot(nus,RY_2)
 subplot(2,2,4)
 plot(nus,FY_2)

tau=taq/2;
expo=exp(-t/tau);
Y_3=Yt.*expo;

sY_3=fftshift(fft(Y_3));
RY_3=real(sY_3);
IY_3=imag(sY_3);
MY_3=abs(sY_3);

figure
subplot(2,2,1)
plot(t,Y)
subplot(2,2,2)
plot(nus,MY_3)
subplot(2,2,3)
plot(t,expo)
subplot(2,2,4)
plot(nus,RY_3)




Amp=[4 2 3 6];
nu = [1000 2500 9000 45000];
phi = [0 pi/2 pi 0]; 

Yt= zeros(1,naq);


for i = 1:length(Amp)
    Yt= Yt+Amp(i)*cos(2*pi*nu(i)*t+phi(i));
end

porte= zeros(1,naq);
porte(2000:8000)=1;

sYt=fftshift(fft(Yt));
fin=porte.*sYt

abs(fin)
plot(nus,abs(fin))

iY=ifft(ifftshift(sYt))

figure
subplot(2,2,1)
plot(t,Yt)
subplot(2,2,2)
plot(nus,abs(sYt))
subplot(2,2,3)
plot(t,iY)
subplot(2,2,4)
plot(nus,abs(fin))



