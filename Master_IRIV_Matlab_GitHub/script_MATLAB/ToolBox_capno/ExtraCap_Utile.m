


%% Selection de la portion Utile de la courbe


function [y,N] = ExtraCap_Utile(y)
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
    N=length(y); % N = longueur du vecteur y après extraction de la portion utile de la courbe
end





