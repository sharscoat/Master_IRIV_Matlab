% % % Courbe type de CAPNOGRAM
% % % 12/06/17, Sebastien HARSCOAT. 




%% initialisation
clc
clear
%% 
close all


% Tfonc='c+((d-c)-a*(t1-t)).*(t>t0).*(t<=t1).*(1-exp(-(t-t0)/T1))+(d-c)*exp(-(t-t1)/T2).*(t>t1)';

% TYfonc='(t>=t0).*(t<t1).*(a/(1-exp(-alpa*t1)))*(1-exp(-alpha*t)) + (t>t1).*(a*(exp(-beta*(t-t1))-exp(-beta*(t2-t1)))/(1-exp(-beta*(t2-t1)))'

% (t>=t0).*(t<t1).*(a/(1-exp(-alpha*t1)))*(1-exp(-alpha*t))+(a*(exp(-beta*(t-t1))-exp(-beta*(t2-t1))))./(1-exp(-beta*(t2-t1))).*(t>t1);


%% essaie de fonction 

t = [0:1:100];

%t0=0;
t1=50;
t2=100;
alpha=0.5;
beta=0.5;
a=0.8;

y=(t<t1).*(a/(1-exp(-alpha*t1))).*(1-exp(-alpha*t))+a*((exp(-beta*(t-t1))-exp(-beta*(t2-t1)))/(1-exp(-beta*(t2-t1)))).*(t>=t1);
% (t>t0).*


figure(1);
plot(t,y,'r','-');




test1_y=(t<t1).*(a/(1-exp(-alpha*t1))).*(1-exp(-alpha*t));
plot(test1_y);

test2_y=a*((exp(-beta*(t-t1))-exp(-beta*(t2-t1)))/(1-exp(-beta*(t2-t1)))).*(t>=t1);
plot(test2_y);


monte=(t<t1).*(a/(1-exp(-alpha*t1))).*(1-exp(-alpha*t));
descente=a*((exp(-beta*(t-t1))-exp(-beta*(t2-t1)))/(1-exp(-beta*(t2-t1)))).*(t>=t1);
figure();
plot(monte+descente);





%%

t = [0:1:100];

% Paramètres
t1=70;
t2=100;
alpha=0.08;
beta=0.5;
a=0.8;

% Equation de la courbe
courbe_1 = (t<t1).*(a/(1-exp(-alpha*t1))).*(1-exp(-alpha*t)) + a*((exp(-beta*(t-t1))-exp(-beta*(t2-t1)))/(1-exp(-beta*(t2-t1)))).*(t>=t1);
monte=(t<t1).*(a/(1-exp(-alpha*t1))).*(1-exp(-alpha*t));
descente=a*((exp(-beta*(t-t1))-exp(-beta*(t2-t1)))/(1-exp(-beta*(t2-t1)))).*(t>=t1);
courbe = monte + descente

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
title(nom_fich);
subplot(3,2,5:6);
plot(1:N,data,'r.-'); 
ylim([0 1.2]); grid on;


%% Figure 3D : Scalogramme en 3D.

figure(4);
title('fichier :');
mesh(resultat,[-5 5]);colorbar;  %





%% Nouvelle fonction de courbe

t = [0:1:100];

% Paramètres
t1=70;
T1=0.18;
T2=0.5;
a=0.8;
d=0;

% équation de la courbe
monte_B = (d-a*(t1-t)).*(1-exp(-t/T1)).*(t<=t1);
descente_B = d*exp(-(t-t1)/T2).*(t>t1);
courbe_B = monte_B + monte_B;

% Figure
figure();
plot(monte_B,'b.-');
