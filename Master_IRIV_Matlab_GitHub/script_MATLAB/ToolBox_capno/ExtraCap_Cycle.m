%% EN fait, on va essayer de trouver les différents cycles ou CAPNOGRAMME

function [B,NP] = ExtraCap_Cycle(y)
    Ymin=min(y);
    Ymax=max(y);
    i=1;
    index=1;
    cmax=0.4;
    cmin=0.1;
    B(index)=i;
    while i<N
        imax =  i     -1 + find( (y(i    :end)-Ymin) > (Ymax-Ymin)*cmax,1, 'first');
        imin1 = imax  -1 + find( (y(imax :end)-Ymin) < (Ymax-Ymin)*cmin,1,'first');
        imin2 = imin1 -1 + find( (y(imin1:end)-Ymin) > (Ymax-Ymin)*cmin,1,'first');
        i=floor((imin1+imin2)/2);
        if isempty(i)
            i=N;
        else
            index=index+1;
            B(index)=i;
            i=imin2;
        end
    end

    NP=length(B)-1;
    
end