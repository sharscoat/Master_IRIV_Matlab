% % % essais capnogrammes avec ondelettes de Haar 
% % % 19/09/17, YT. 




%% initialisation
clc
clear
%% 
close all



%% sélection du fichier

[F,Par]=uigetfile('*.txt','FICHIER A OUVRIR','/Users/sharscoat/Documents/SHARSCOAT/Master_IRIV/Capnographie/CAPNO_ref_propre');
filename=[Par F];

% findstr
% strfind(vect,':')
% NOM du FICHIER
n_fich=filename(end-15:end);
debut=strfind(n_fich,'/')+1;
fin=strfind(n_fich,'.')-1;
nom_fich=n_fich(debut:fin);


%% ouverture, conversion et récupération des données

FILE=fopen(filename);
T=textscan(FILE,'%s');
fclose(FILE);
T=T{1};
T=strrep(T,',','.');
T(1:10)=[];
T(end-2:end)=[];
c=T(2:2:end);
y=cellfun(@str2num,(c(1:end)));
tps=[0:1:length(y)-1];
t=tps;


%% Selection de la portion Utile de la courbe

Y=max(y);

N=length(y);
i1=find(y>Y/20,1,'first')-10;
i2=find(y>Y/20,1,'last')+10;
if i1<1
    i1=1;
end
if i2>N
    i2=N;
end

if i1>1
    t(1:(i1-1))=[];
    y(1:(i1-1))=[];
end
if i2<N
    t(i2+1:end)=[];
    y(i2+1:end)=[];
end
N=length(y);

%%% la  portion utile de la courbe est désormais choisie




%% ANALYSE SCALOGRAMME

% donnees 
data=y; 
data=data/max(abs(data)); 
N=length(data); 
figure(1); 
plot(1:N,data,'r.-'); 
ylim([0 1.5]); grid on; 

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
ylim([0 1.5]); grid on;








