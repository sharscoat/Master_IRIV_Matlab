%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sébastien HARSCOAT le 20/06/2017
% Programme : "courbe_CAPNO"
% permet de produire une courbes de CAPNOGRAM et de réaliser son SCALOGRAM
% associé.



%% Fonction : "courbe_CAPNO" 


function[courbe,resultat] = courbe_CAPNO(t1,t2, T1, T2, a, d)

    t = [0:1:t2];    % t = [0:1:100];
    t0=0;
    c=0;

    % Paramètres
    % t0=0; t1=70; t2=100;
    % T1=10;   % les valeur de "T1" son comprise entre 2 et 20
    % T2=5;
    % a=0.003;  % les valeur de "a" son comprise entre 0 et 0.1
    % c=0; d=5;

    % Equation de la courbe
    courbe_1 = c + ((d-c)-a*(t1-t)).*(t>t0).*(t<=t1).*(1-exp(-(t-t0)/T1)) + (d-c)*exp(-(t-t1)/T2).*(t>t1);
    monte =((d-c)-a*(t1-t)).*(t>t0).*(t<=t1).*(1-exp(-(t-t0)/T1));
    descente =(d-c)*exp(-(t-t1)/T2).*(t>t1);
    courbe = monte + descente;

    
    % SCALOGRAM
    data=courbe';
    data=data/max(abs(data)); 
    N=length(data); 
        %figure(1); 
        %plot(1:N,data,'r.-'); 
        %ylim([0 1.2]); grid on; 

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

end














