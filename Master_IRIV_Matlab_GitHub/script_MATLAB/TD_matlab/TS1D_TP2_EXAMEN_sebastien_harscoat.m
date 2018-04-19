clear all
close all


% =========================================================================================
% PARTIE 1
% Paramètres d'échantillonnage du signal et création du signal
nue = 1.0000e+05;     % fréquence d'échantillonnage
te = 1/nue;      % période d'échantillonnage
nacq = 10001;	% nombre d'échantillons
tacq = (nacq-1)*te;	% durée d'acquisition
t = [0:te:tacq];       % vecteur des temps d'acquisition
dnu = 1/tacq;     % résolution fréquentielle
numax = nue/2;   % fréquence maximale
nus = [-numax:dnu:numax];     % vecteur des fréquences

% écrire une fonction f(t) exemple somme de 3 cosinus
%%% Somme de fonctions sinusoidales
%             Amp1 = 3;
%             nu1 = 1000;
%             phi1 = pi/2;
%             Y1 = Amp1*cos(2*pi*nu1*t + phi1);
% 
%             Amp2 = 2;
%             nu2 = 4000;
%             phi2 = 0;
%             Y2 = Amp2*cos(2*pi*nu2*t + phi2);
% 
%             Amp3 = 4;
%             nu3 = 8000;
%             phi3 = pi;
%             Y3 = Amp3*cos(2*pi*nu3*t + phi3);
% 
%             y = Y1 + Y2 + Y3;
            % ou encor

y = 3*cos(2*pi*1000*t + pi/2) + 2*cos(2*pi*4000*t + 0) + 4*cos(2*pi*8000*t + pi);

% calcul du spectre
sy = fftshift(fft(y));

% partie réel
Ry = real(sy);
% partie imaginaire
Iy = imag(sy);
% module 
My = abs(sy);

% représentation graphique
figure(1);
subplot(2,2,1)
plot(t(1:500),y(1:500)); % signal
    title('Signal');
    set(get(gca,'XLabel'),'String','Temps (s)','FontName','Arial','FontSize',9)
    set(get(gca,'YLabel'),'String','Amplitude (V)','FontName','Arial','FontSize',9)
subplot(2,2,2)
plot(nus,Ry);	% partie réelle du spectre
    title('partie réelle');
    set(get(gca,'XLabel'),'String','Temps (s)','FontName','Arial','FontSize',9)
    set(get(gca,'YLabel'),'String','Fréquence (Hz)','FontName','Arial','FontSize',9)
subplot(2,2,3)
plot(nus,Iy);	% partie imaginaire du spectre
    title('partie imaginaire');
    set(get(gca,'XLabel'),'String','Temps (s)','FontName','Arial','FontSize',9)
    set(get(gca,'YLabel'),'String','Fréquence (Hz)','FontName','Arial','FontSize',9)
subplot(2,2,4)
plot(nus,My);	% module du spectre
    title('partie module');
    set(get(gca,'XLabel'),'String','Temps (s)','FontName','Arial','FontSize',9)
    set(get(gca,'YLabel'),'String','Fréquence (Hz)','FontName','Arial','FontSize',9)
   

% ========================================================================================
% PARTIE 2
% écriture d'une fonction porte et multiplication du cosinus par porte (temps de mesure) 

% définition de la fonction porte
porte = zeros(1,nacq);
porte(1:(nacq-1)/20) = 1;


% Filtrage du signal
yp= y.*porte ;

% calcul du spectre
syp = fftshift(fft(yp));

% partie réel
Ryp = real(syp);
% partie imaginaire
Iyp = imag(syp);
% module 
Myp = abs(syp);

% représentation graphique
figure(2);
subplot(2,2,1)
plot(t(1:1000),yp(1:1000)); % signal filtré
subplot(2,2,2)
plot(nus,Ryp);	% partie réelle du spectre
subplot(2,2,3)
plot(nus,Iyp);	% partie imaginaire du spectre 
subplot(2,2,4)
plot(nus,Myp);	% module du spectre

% ========================================================================================
% PARTIE 3
% filtrage par une fonction améliorée

%définition de la fonction de filtrage
tau = tacq/20;
filtre = exp(-t/tau);

% Filtrage du signal
yf= y .* filtre;

% calcul du spectre
syf = fftshift(fft(yf));

% partie réel
Ryf = real(syf);
% partie imaginaire
Iyf = imag(syp);
% module 
Myf = abs(syf);


% représentation graphique
figure(2);
subplot(2,2,1)
plot(t(1:5000),yf(1:5000)); % signal filtré
subplot(2,2,2)
plot(nus,Ryf);	% partie réelle du spectre
subplot(2,2,3)
plot(nus,Iyf);	% partie imaginaire du spectre 
subplot(2,2,4)
plot(nus,Myf);	% module du spectre


% ========================================================================================
% PARTIE 4
% génération du signal bruité

% ajout du bruit
namp = 5;    % amplitude du bruit
sn = namp * rand(nacq,1);  % bruit
yn= y .* sn';   % signal bruité

%calcul du spectre
syn= fftshift(fft(yn));

% module 
Myn = abs(syn);

% représentation graphique
figure(4);
subplot(2,2,1); 
plot(t(1:250),y(1:250))   % signal non bruité
subplot(2,2,2); 
plot(t(1:250),yn(1:250)); % signal bruité
subplot(2,2,3)
plot(nus,My);  % module du spectre du signal non bruité
subplot(2,2,4)
plot(nus,Myn);   % module du spectre du signal bruité



% =========================================================================================
% PARTIE 5
% Filtre RC simple

% écrire la constante de temps RC avec une fréquence de coupure à 
rc = 5*te/(2*pi); %fréquence de coupure à 4000Hz si nue = 20000Hz 
nue = 20000;
te = 1/nue;
nuc = nue/5


% écrire la réponse impulsionnelle
repimp= exp(-t/rc)/rc;

% écrire la fonction de transfer 
trans= 1/(1 + 2*pi*j*nue*rc);

% écrire gain du filtre
gain= -log10(1 + (2*pi*nuc*rc)^2);

%écrire le déphasage du filtre
phase= arctan(-2*pi*nue*rc) ;  

% TF de la fonction de transfert h(t) --> TF --> H(nu)
Hnu = fftshift(fft(trans));

% représentation graphique
figure(5);
subplot(2,2,1)
plot(t,repimp); % Réponse impulsionnelle
axis([-0.05 0.1 0 3000]);
subplot(2,2,2)
plot(nus,trans); % Fonction de transfert
subplot(2,2,3)
plot(nus,gain); % title('Gain (dB)')
subplot(2,2,4)
plot(nu,Hnu); % Fonction de transfert


% ==========================================================================================
% PARTIE 6
% Filtre RC simple: Filtrage temporel

% écrire la convolution du signal bruité avec la réponse impulsionnelle
yc = conv(yn,repimp); %utilisez la fonction conv de Matlab
yc2=yc(1:nacq); %adaptation de l'échelle apres convolution

% calcul du spectre
syc2 = fftshift(fft(yc2));

% module 
Myc2 = abs(syc2);

% représentation graphique
figure(6);
subplot(2,2,1)
plot(t(1:250),yn(1:250)); %signal bruité
subplot(2,2,2)
plot(t,yc2) %signal convolué
subplot(2,2,3)
plot(nus,Myn); %spectre signal non filtré
subplot(2,2,4)
plot(nus,Myc2); %spectre signal filtré par convolution


% ==========================================================================================
% PARTIE 7
% Filtre RC simple: Filtrage fréquentiel

% écrire le signal filtré dans le domaine fréquentiel
syc3 = fftshift(fft(yn) * trans);

% module 
Myc3 = abs(syc3);


% représentation graphique
figure(7);
subplot(2,2,1)
plot(t(1:250),yn(1:250)); %signal bruité
subplot(2,2,2)
plot(t,yfc) %spectre signal non filtré
subplot(2,2,3)
plot(); %spectre signal filtré par convolution
subplot(2,2,4)
plot(); %spectre signal filtré dans domaine fréquentiel







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%    OPTIONNEL
%===========================================================================================
% Fonctions Matlab évoluées (Signal Processing Toolbox)
% Entrée manuelles des paramètres du filtre
% titre= 'Select filter';
% params = {'ordre','f coupure'};
% lines =2;
% defaut= {'6','0.5'};
% entrees = inputdlg(params,titre,lines,defaut);
% [ns]= deal(entrees{1});
% n = str2num(ns);
% [cs]= deal(entrees{2});
% co = pi*str2num(cs)/5;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construction d'un filtre Butterworth
n = 6; % ordre du filtre
co = .4; % fréquence de coupure relative(nuc/nushannon)
[num1,den1] = butter(n,co); % Butterworth
%filtrage du signal
yfb=filter(num1,den1,yn);
%Calculez le spectre
syfb= ;



% représentation graphique
figure(8);
subplot(2,2,1)
plot() %signal bruité
subplot(2,2,2)
plot(); %signal filtré
subplot(2,2,3)
plot(); %spectre signal non filtré
subplot(2,2,4)
plot(); %spectre signal filtré


%%%% Autres filtres à essayer:
%[num1,den1] = cheby1(n,.5,co); % Chebychev
%[num1,den1] = butter(n,[co/8 co],'bandpass');



