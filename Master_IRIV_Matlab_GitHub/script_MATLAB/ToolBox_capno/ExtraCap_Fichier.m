
%  ouverture, conversion et récupération des données

function [y,t] = ExtraCap_Fichier(filename)
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
end 
