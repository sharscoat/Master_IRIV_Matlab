%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%    Essai de Programme pour créer une data base    %%%%%%%%%%%%%
%%%%%%%%%%% PROGRAMME ANALYSE COURBE CAPNOGRAPHIE - SCALOGRAM %%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ''' Essaie pour basculer des données de type "data base" dans un fichier
% utilisation par "R". Par la suite analyse des différentes courbes "CO2",
% "Flow" et "O2" données dans les fichiers.


%% initialisation
clc
clear
close all


%% récuperer les noms des fichiers dans un dossier

%[F,Par]=uigetfile('*.txt','FICHIER A OUVRIR');
%filename=[Par F];

file_dossier='/Users/sharscoat/Documents/SHARSCOAT/Master_IRIV/Capnographie/0_problemes/209786.txt'
[F,Par]=uigetfile(file_dossier,'FICHIER A OUVRIR');
filename=[Par F];


%% ouverture, conversion et récupération des données
FILE=fopen(filename);
T=textscan(FILE,'%s');
fclose(FILE);
T=T{1};
T=strrep(T,',',',');
T(1:12)=[];
T(end-2:end)=[];
% pour le "CO2"
c_co2=T(2:4:end);
y_co2=cellfun(@str2num,(c_co2(1:end)));
% pour le "Flow"
c_flow=T(3:4:end);
y_flow=cellfun(@str2num,(c_flow(1:end)));
% pour le "Flow"
c_o2=T(4:4:end);
y_o2=cellfun(@str2num,(c_o2(1:end)));

% pour y
c=T(2:4:end);
y=cellfun(@str2num,(c(1:end)));

% pour t
tps=[0:1:length(y_co2)-1];
t=tps;



%% Créer une DATA_BASE avec les différentes variables : "y_co2", "y_flow" et "y_O2"

data_base = table(y_co2, y_flow, y_o2);    % création d'une base de données comme dans "R"

data_base(1:20,:);       % visualisation des 20 premières lignes de la base de données "data_base"

data_base.y_co2(1:20,:);  % visualisation des 20 premières lignes de la variable "y_co2" de la base "data_base"


%% enregistrement des variables "y_co2", "y_flow" et "y_O2" dans un fichier ".m"

% ne permet que d'enregistrer une variable ou un vecteur ou une matrice
save('/Users/sharscoat/Documents/SHARSCOAT/Master_IRIV/Capnographie/essaie_R.m','y_co2','y_flow','y_o2')

save('/Users/sharscoat/Documents/SHARSCOAT/Master_IRIV/Capnographie/essaie_R.txt','y_co2','y_flow','y_o2', '-ascii')


% enregsitrement de la base de données "data_base" dans un fichier ".txt"
% permet d'enregistrer ou plutôt d'écrire une base de données ici
% "data_base" dans un fichier ".txt"

writetable(data_base,'/Users/sharscoat/Documents/SHARSCOAT/Master_IRIV/Capnographie/essaie_database_R.txt');



%% utilisation du fichier ".txt" dans "R"

% data_base_CO2 <- read.table("/Users/sharscoat/Documents/SHARSCOAT/Master_IRIV/Capnographie/essaie_database_R.txt",header = T,sep = ',',dec = ".")










%% Graphique 

figure();

subplot(3,1,1);
plot(data_base.y_co2,'b.-');
grid on;

subplot(3,1,2);
plot(data_base.y_flow, '-', 'Color', [0.87 0.49 0]);
grid on;

subplot(3,1,3);
plot(data_base.y_o2,'r.-');
grid on;



%%%%%%%%%%%%%%%%%%%%%%
figure();

subplot(3,1,1);
plot(data_base.y_co2(1000:end),'b.-');
grid on;

subplot(3,1,2);
plot(data_base.y_flow(1000:end), '-', 'Color', [0.87 0.49 0]);
grid on;

subplot(3,1,3);
plot(data_base.y_o2(1000:end),'r.-');
grid on;



%% Dérivée Première

tps=[0:1:length(y_co2)-1];
t=tps';

% pour "y_co2"
dy_co2 = diff(data_base.y_co2)/diff(t);
plot(dy_co2);

% pour "y_flow"
dy_flow = diff(data_base.y_flow)/diff(t);
plot(dy_flow);

% pour "y_co2"
dy_o2 = diff(data_base.y_o2)/diff(t);
plot(dy_o2);



% toutes les dérivées premières
figure();

subplot(3,1,1);
plot(dy_co2,'b.-');
grid on;

subplot(3,1,2);
plot(dy_flow, '-', 'Color', [0.87 0.49 0]);
grid on;

subplot(3,1,3);
plot(dy_o2,'r.-');
grid on;




% graphique courbe "Flow" ou "Débit" avec sa Dérivées Premières

figure();

subplot(2,1,1);
plot(data_base.y_flow, '-', 'Color', [0.87 0.49 0]);
grid on;

subplot(2,1,2);
plot(dy_flow, '-');
grid on;



% graphique courbe "O2" avec sa Dérivées Premières

figure();

subplot(2,1,1);
plot(data_base.y_o2,'r.-');
grid on;

subplot(2,1,2);
plot(dy_o2, '-');
grid on;








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

    %%% la  portion utile de la courbe est désormais choisie



    %% EN fait, on va essayer de trouver les différents cycles sans findpeaks

    Ymin=min(y);
    Ymax=max(y);
    i=1;
    index=1;
    cmax=0.4;
    cmin=0.1;
    B = [];
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

    NP = 0;
    NP=length(B)-1;
    
    if B(end)>length(y)
        B(end)=length(y);
    end
    disp(B')
    disp(NP)

    %% indexation terminée
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    
 %% Graphique 

figure();

subplot(3,1,1);
plot(data_base.y_co2,'b.-');
hold on
plot(t(B),y(B),'sk','MarkerFaceColor','k', 'MarkerSize',8)
grid on;

subplot(3,1,2);
plot(data_base.y_flow, '-', 'Color', [0.87 0.49 0]);
hold on
plot(t(B),data_base.y_flow(B),'sk','MarkerFaceColor','k', 'MarkerSize',8)
grid on;

subplot(3,1,3);
plot(data_base.y_o2,'r.-');
hold on
plot(t(B),data_base.y_o2(B),'sk','MarkerFaceColor','k', 'MarkerSize',8)
grid on;
    
    
    
    
    
    
    
