%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% S�bastien HARSCOAT le 16/06/2017
% Programme : "load_CAPNO"
% permet de charger des courbes de capnographie avec le nom de fichier



%% Fonction : "load_CAPNO"


function[y] = load_CAPNO()
    % s�lection du fichier
    [F,Par]=uigetfile('*.txt','FICHIER A OUVRIR','/Users/sharscoat/Documents/SHARSCOAT/Master_IRIV/Capnographie/CAPNO_ref_propre');
    filename=[Par F];

    n_fich=filename(end-15:end);
    debut=strfind(n_fich,'/')+1;
    fin=strfind(n_fich,'.')-1;
    nom_fich=n_fich(debut:fin);

    % ouverture, conversion et r�cup�ration des donn�es
    FILE=fopen(filename);
    T=textscan(FILE,'%s');
    fclose(FILE);
    T=T{1};
    T=strrep(T,',','.');
    T(1:10)=[];
    T(end-2:end)=[];
    c=T(2:2:end);
    y=cellfun(@str2num,(c(1:end)));
    
end





