%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Courbe de CAPNOGRAM et Son SCALOGRAM
% Sébastien HARSCOAT le 13/06/2017



%% initialisation
clc
clear
%% 
close all



%%  Construction de courbe de type "EXPONENTIELLE"

t = [0:1:100];

% Paramètres
t1=70;
t2=100;
alpha=0.70;
beta=0.5;
a=0.8;

% Equation de la courbe
courbe_1 = (t<t1).*(a/(1-exp(-alpha*t1))).*(1-exp(-alpha*t)) + a*((exp(-beta*(t-t1))-exp(-beta*(t2-t1)))/(1-exp(-beta*(t2-t1)))).*(t>=t1);
monte=(t<t1).*(a/(1-exp(-alpha*t1))).*(1-exp(-alpha*t));
descente=a*((exp(-beta*(t-t1))-exp(-beta*(t2-t1)))/(1-exp(-beta*(t2-t1)))).*(t>=t1);
courbe = monte + descente;

% Figure
figure();
plot(courbe_1,'b.-');


% SCALOGRAM
data=courbe_1';
data=data/max(abs(data)); 
N=length(data); 
figure(1); 
plot(1:N,data,'r.-'); 
ylim([0 1.2]); grid on; 

% conversion en puissance de 2 
n=floor(log(N)/log(2)); 
N0=2^n; 

% calcul de scalogrammes 
mask_l=[1 1]/sqrt(2); 
mask_h=[1 -1]/sqrt(2); 
scalogramme=data(1:N0); 
buffer=data(1:N0); 
for k=n:-1:1, 
    % % mere 
    data_h=filter(mask_h,1,buffer); 
    data_k=data_h(1:2:end); 
    % % pere 
    data_l=filter(mask_l,1,buffer); 
    buffer=data_l(1:2:end); 
    
    % % resizing 
    temp=kron(data_k,ones(2^(n-k+1),1)); 
    scalogramme=[scalogramme temp]; 
end 
resultat=[fliplr(scalogramme)]'; 


% affichage 
figure(2); 
imagesc(resultat,[-5 5]); colorbar; 
xlabel('n (temps)'); 
ylabel('k (echelle)'); 



%% Figure : Scalogramme et Courbe de Capnogramme




figure(3);
title('fichier :');
subplot(3,2,1:4);
imagesc(resultat,[-5 5]); colorbar; 
xlabel('n (temps)'); 
ylabel('k (echelle)');
subplot(3,2,5:6);
plot(1:N,data,'r.-'); 
ylim([0 1.2]); grid on;


%% Figure 3D : Scalogramme en 3D.

figure(4);
title('fichier :');
mesh(resultat,[-5 5]);colorbar;  %



%% Figure résumé

figure(5);
title('fichier :');
subplot(2,2,1);
imagesc(resultat,[-5 5]); colorbar; 
xlabel('n (temps)'); 
ylabel('k (echelle)');

subplot(2,2,3);
plot(1:N,data,'r.-'); 
ylim([0 1.2]); grid on;

subplot(2,2,2);
mesh(resultat,[-5 5]);colorbar; 


imagesc(A,[-5 5]); caxis([0,0.5]); colorbar;

imagesc(A,[-5 5]); caxis([0,2]); colorbar;imagesc(A,[-5 5]); caxis([0,2]); colormap('jet');colorbar;






%%  Construction de courbe de type "EXPONENTIELLE + DROITE"

% 'c+((d-c)-a*(t1-t)).*(t>t0).*(t<=t1).*(1-exp(-(t-t0)/T1)) + (d-c)*exp(-(t-t1)/T2).*(t>t1)';

t = [0:1:100];

% Paramètres
t0=0
t1=70;
t2=100;
T1=10;   % les valeur de "T1" son comprise entre 2 et 20
T2=5;
a=0.003;  % les valeur de "a" son comprise entre 0 et 0.1
c=0;
d=5;

% Equation de la courbe
courbe_1 = c + ((d-c)-a*(t1-t)).*(t>t0).*(t<=t1).*(1-exp(-(t-t0)/T1)) + (d-c)*exp(-(t-t1)/T2).*(t>t1);
monte=((d-c)-a*(t1-t)).*(t>t0).*(t<=t1).*(1-exp(-(t-t0)/T1));
descente=(d-c)*exp(-(t-t1)/T2).*(t>t1);
courbe = monte + descente;

% Figure
figure();
plot(courbe_1,'b.-');
plot(courbe,'b.-'); ylim([0 6]); grid on;

