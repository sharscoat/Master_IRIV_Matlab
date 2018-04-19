%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%              TD_1   Traitement de Signal 1D                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
close all

%=========================================================================
% Partie 1

%%% vecteur temps
% pr�cisions te
te = 0.00001;
% nombre acquisition
nacq = 10001;
% borne final acquisition
tacq = (nacq-1)*te;
% vecteur temps t
t = [0:te:tacq];



%%% Cr�ation d'un signal 

Amp = 1000;
nu = 1000;
phi = 0;

% fonction du signal : Y = Amp*cos(2*pi*nu*t + phi)
Y = Amp*cos(2*pi*nu*t + phi);

figure;
plot(t,Y);



%====================================================

%%% Transform�e de Fourier de notre signal
sY = fft(Y);
plot(sY);
sY = fftshift(fft(Y));
plot(sY);

% partie r�el
RY = real(sY);
% partie imaginaire
IY = imag(sY);
% module 
MY = abs(sY);

plot(RY);

% fr�quence d'�chantillonnage 'nue'
nue = 1/te;
%
numax = nue/2;
% incr�mentation 'dnu'
dnu = 1/tacq;
% vecteur de fr�quence
nus = [-numax:dnu:numax];

plot(nus,MY);





%====================================================

%%% Somme de fonctions sinusoidales
Amp1 = 4;
nu1 = 1000;
phi1 = 0;
Y1 = Amp1*cos(2*pi*nu1*t + phi1);

Amp2 = 2;
nu2 = 2500;
phi2 = pi/2;
Y2 = Amp2*cos(2*pi*nu2*t + phi2);

Amp3 = 3;
nu3 = 9000;
phi3 = pi;
Y3 = Amp3*cos(2*pi*nu3*t + phi3);

Y = Y1 + Y2 + Y3;


%%% Transform�e de Fourier
sY = fft(Y);
sY = fftshift(fft(Y));
plot(sY);

% partie r�el
RY = real(sY);
% partie imaginaire
IY = imag(sY);
% module 
MY = abs(sY);

plot(RY);

% fr�quence d'�chantillonnage 'nue'
te = 0.00001;
nue = 1/te;
%
numax = nue/2;
% incr�mentation 'dnu'
dnu = 1/tacq;
% vecteur de fr�quence
nus = [-numax:dnu:numax];

plot(nus,MY);


%%%%%%%%%%%%   avec une fonction %%%%%%%%%%%%%
Amp = [4 2 3];
nu = [1000 2500 9000];
phi = [0 pi/2 pi];
Y = zeros(1,nacq);

for i = 1:length(Amp)
    Y = Y + Amp(i)*cos(2*pi*nu(i)*t + phi(i));
end


% partie r�el
RY = real(sY);
% partie imaginaire
IY = imag(sY);
% module 
MY = abs(sY);
% phase
AY = unwrap(angle(sY));   % 'unwrap' permet de contraidre entre 0 et 2*pi

figure;
subplot(2,2,1);
plot(t(1:300),Y(1:300));
subplot(2,2,2);
plot(nus,MY);
subplot(2,2,3);
plot(nus,RY);
subplot(2,2,4);
plot(nus,AY);




%%%% Fr�quence de Nyquist
% si on choisi mal les fr�quences d'�chantillonnage on risque d'avoir un
% replie spectral.
% on doit prendre des fr�quence d'�chantillonnage 2 fois sup�rieur au
% minimum de la fr�quence de notre signal.







%====================================================

% on cherche � multiplier une fonction porte de valeur 1 entre 0 et
% (nacq-1)/20 et valant 0 pour le reste des valeurs.

%%% cr�ation d'une fonction porte :

nacq = 10001
porte = zeros(1,nacq);

for i = 1:(nacq-1)/20
    porte(i) = 1;
end

for i = ((nacq-1)/20)+1:nacq
    porte(i) = 0;
end

% une autre fa�on de faire
nacq = 10001
Te = 60
porte = zeros(1,nacq);
porte(1:(nacq-1)/Te) = 1;


%%% multiplication de la fonction 'porte' par la fonction Y.
Yt = Y.*porte;   % c'est un produit point � point (et non un produit matriciel) donc .*


%%% Transform�e de Fourier

% fr�quence d'�chantillonnage 'nue'
te = 0.00001;
nue = 1/te;
%
numax = nue/2;
% incr�mentation 'dnu'
dnu = 1/tacq;
% vecteur de fr�quence
nus = [-numax:dnu:numax];

sYt = fftshift(fft(Yt));

% partie r�el
RYt = real(sYt);
% partie imaginaire
IYt = imag(sYt);
% module 
MYt = abs(sYt);
% phase
AYt = unwrap(angle(sYt));

figure;
subplot(2,2,1);
plot(t,Yt);
subplot(2,2,2);
plot(nus,MYt);
subplot(2,2,3);
plot(nus,RYt);
subplot(2,2,4);
plot(nus,AYt);


figure;
subplot(2,2,1);
plot(t,Y);
subplot(2,2,2);
plot(nus,MY);
subplot(2,2,3);
plot(t,Yt);
subplot(2,2,4);
plot(nus,MYt);




%====================================================

%%% cr�ation d'une fonction Exponentielle D�croissante : expo = exp(-t/tau)
tau = tacq/60;
expo = exp(-t/tau);


Ye = Y .* expo;   % c'est un produit point � point donc .*



%%% Transform�e de Fourier

% fr�quence d'�chantillonnage 'nue'
te = 0.00001;
nue = 1/te;
%
numax = nue/2;
% incr�mentation 'dnu'
dnu = 1/tacq;
% vecteur de fr�quence
nus = [-numax:dnu:numax];

sYe = fftshift(fft(Ye));

% module 
MYe = abs(sYe);

figure;
subplot(3,2,1);
plot(t(1:1000),Y(1:1000));
subplot(3,2,2);
plot(nus,MY);
subplot(3,2,3);
plot(t(1:1000),Ye(1:1000));
subplot(3,2,4);
plot(nus,MYe);
subplot(3,2,5);
plot(t(1:1000),Yt(1:1000));
subplot(3,2,6);
plot(nus,MYt);






%====================================================

%%% cr�ation d'un filtre

nacq = 10001;           % nombre d'acquisition
te = 0.00001;           % temps entre chaque acquisition
tacq = (nacq-1)*te;     % temps d'acquisition
t = [0:te:tacq];        % 


Amp = [4 2 3 6];
nu = [1000 2500 9000 45000];
phi = [0 pi/2 pi 0];
Y = zeros(1,nacq);

for i = 1:length(Amp)
    Y = Y + Amp(i)*cos(2*pi*nu(i)*t + phi(i));
end


%%% fonction porte dans le domaine des fr�quences 



%%% Transform�e de Fourier

% fr�quence d'�chantillonnage 'nue'
te = 0.00001;
nue = 1/te;
numax = nue/2;
dnu = 1/tacq;
nus = [-numax:dnu:numax];

sY = fftshift(fft(Y));


%%% fonction porte dans le domaine des fr�quences
filtre = zeros(1,nacq);
filtre(2000:8000) = 1;

% on filtre
sYf = sY .* filtre   % on multiplie la fonction 'filtre' dans le domaine de Fourier
                     % donc sur la transform�e de Fourier.



% module 
MY = abs(sY);
MYf = abs(sYf);

figure;
subplot(2,1,1);
plot(nus,MY);
subplot(2,1,2);
plot(nus,MYf);


%%% Transform�e de Fourier inverse  (pour retrouver la fonction d�bruit�e
%%% dans le domaine temporel).

Yf = ifft(ifftshift(sYf));



figure;
subplot(2,2,1);
plot(t,Y);
subplot(2,2,2);
plot(nus,MY);
subplot(2,2,3);
plot(t,Yf);
subplot(2,2,4);
plot(nus,MYf);






%====================================================


