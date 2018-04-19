%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROGRAMME Courbe_CAPNO et SCALOGRAM
% différence entre Scalogram
% Courbe de CAPNOGRAM et Son SCALOGRAM
% Sébastien HARSCOAT le 15/06/2017



%% initialisation
%clc
%clear
%close all



%% Programme "courbe_CAPNO" 

function[courbe] = courbe_CAPNO(t1,t2, alpha, beta, a)

    t = [0:1:t2];    % t = [0:1:100];

    % Paramètres
    % t1=70; t2=100; alpha=0.70; beta=0.5; a=0.8;

    % Equation de la courbe
    monte=(t<t1).*(a/(1-exp(-alpha*t1))).*(1-exp(-alpha*t));
    descente=a*((exp(-beta*(t-t1))-exp(-beta*(t2-t1)))/(1-exp(-beta*(t2-t1)))).*(t>=t1);
    [courbe] = monte + descente;

end



%% Programme "scalogram_CAPNO

function [resultat] = scalogram_CAPNO(data)
    
    % SCALOGRAM
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
    [resultat]=[fliplr(scalogramme)]';

end





%% 

%% Programme "global_CAPNO" 

function[courbe,resultat] = global_CAPNO(t1,t2, alpha, beta, a)

    t = [0:1:t2];    % t = [0:1:100];

    % Paramètres
    % t1=70; t2=100; alpha=0.70; beta=0.5; a=0.8;

    % Equation de la courbe
    monte=(t<t1).*(a/(1-exp(-alpha*t1))).*(1-exp(-alpha*t));
    descente=a*((exp(-beta*(t-t1))-exp(-beta*(t2-t1)))/(1-exp(-beta*(t2-t1)))).*(t>=t1);
    courbe = monte + descente;
    [courbe] = courbe

    
    % SCALOGRAM
    data=courbe';
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
    [resultat]=[fliplr(scalogramme)]';

end










