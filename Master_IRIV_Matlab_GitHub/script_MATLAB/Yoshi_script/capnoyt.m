% % % essais capnogrammes avec ondelettes de Haar 
% % % 19/09/17, YT. 

% init 
% clear all; close all; 

% donnees 
data=load('WaveForm_NORMAL_25.txt'); 
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


