%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sébastien HARSCOAT le 16/06/2017
% Programme : "scalogram_CAPNO"
% permet de réaliser le scalogram




%% Fonction : "scalogram_CAPNO"


function[resultat] = scalogram_CAPNO(data)

% SCALOGRAM
    data=data';
    data=data/max(abs(data)); 
    N=length(data); 
    
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
    %resultat_1=[fliplr(scalogramme)]';
    %[resultat] = resultat_1/sqrt(N0);

end 


