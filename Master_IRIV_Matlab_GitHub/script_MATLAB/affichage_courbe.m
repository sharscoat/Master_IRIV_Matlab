%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                                                               %%%%
%%%%                       Script de Travail                       %%%%
%%%%               CAPNOGRAPHY - affichage des courbes             %%%%
%%%%                                                               %%%%
%%%%                       Sébastien Harscoat                      %%%%
%%%%                                                               %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% initialisation
clc
clear
%% 
close all



%% sélection du fichier

[F,Par]=uigetfile('*.txt','FICHIER A OUVRIR');
filename=[Par F];
% filename2=[ filename(1:end-4) '_.csv'];
% filename3=[ filename(1:end-4) '.pdf'];
% filename4=[ filename(1:end-4) '_Res.csv'];
% '/Users/sharscoat/Documents/SHARSCOAT/Master_IRIV/Capnographie/MatLab/scripts/10581GUHMANN_EMPHYSEME.txt'

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






%% ENREGISTREMENT DE LA FIGURE  (sans afficher le graphique)

h=figure();(set(h, 'Visible', 'off') );   % évite d'afficher l'image !!!
plot(y);

saveas(h,'/Users/sharscoat/Documents/SHARSCOAT/Master_IRIV/Capnographie/image_y.jpg','jpg');


%imwrite(plot(y),'/Users/sharscoat/Documents/SHARSCOAT/Master_IRIV/Capnographie/image_y.jpg','jpg','Quality',100);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%        CELA FONCTIONNE !!!!!!!        %%%%%%%


%% récuperer les noms des fichiers dans un dossier

%[F,Par]=uigetfile('*.txt','FICHIER A OUVRIR');
%filename=[Par F];

file_dossier='/Users/sharscoat/Documents/SHARSCOAT/Master_IRIV/Capnographie/0_problemes/209786.txt'
[F,Par]=uigetfile(file_dossier,'FICHIER A OUVRIR');


path = Par ; % répertoire courant
D = dir(path); % récupère un tableau de structure
D = D(~cell2mat({D(:).isdir})); % filter pour ne garder que les noms de fichiers
Liste = {D(:).name}; % transformer en un tableau de cellules texte
Liste = Liste';
Liste_fich = Liste(2:length(Liste));

% tous les chemin pour chaque fichier :

for i = 1:length(Liste_fich)
    liste_path{i}=[Par Liste_fich{i}];
end


% "Liste_fich" donne un liste de nom de fichier se finissant par ".txt"
% Liste_fich{j}(1:end-4); permet de récupérer le nom du fichier sans ".txt"




% tous les nom de fichier avec ".jpg" pour chaque fichier :

for i = 1:length(Liste_fich)
    liste_jpg{i}=[Liste_fich{i}(1:end-4) '.jpg'];
end




%% Enregistrement des graphiques dans un dossier

[F2,Par2]=uigetfile('*.txt','FICHIER A OUVRIR');
filename_res2=[Par2 F2];


% tous les chemin pour enregister les fichier :

for i = 1:length(liste_jpg)
    liste_path2{i}=[Par2 liste_jpg{i}];
end



%% Analyse de tous les fichiers

for j = 1:length(liste_path)
    %% ouverture, conversion et récupération des données
    %filename = liste_path{j};
    dir(liste_path{j})
    FILE=fopen(liste_path{j});
    T=textscan(FILE,'%s');
    fclose(FILE);
    T=T{1};
    T=strrep(T,',',',');
    T(1:12)=[];
    T(end-2:end)=[];
    c=T(2:4:end);
    y=cellfun(@str2num,(c(1:end)));
    tps=[0:1:length(y)-1];
    t=tps;
    


    %% Selection de la portion Utile de la courbe

    ymin = min(y);
    if ymin>0
        y = y - ymin;
    else
        y = y + abs(ymin);
    end
    
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
    
    %% ENREGISTREMENT DE LA FIGURE  (sans afficher le graphique)
    h=figure();(set(h, 'Visible', 'off') );   % évite d'afficher l'image !!!
    plot(y);

    saveas(h,liste_path2{j},'jpg');
    
end












