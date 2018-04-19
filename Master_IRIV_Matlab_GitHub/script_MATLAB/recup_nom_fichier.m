%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% récuperer les noms des fichiers dans un dossier


clc
clear
close all

[F,Par]=uigetfile('*.txt','FICHIER A OUVRIR');
filename=[Par F];

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



for i = 1:length(Liste_fich)
    disp([Par Liste_fich{i}]);
end



%%

[Par,Liste{2}];

% tous les chemin pour chaque fichier :

 
for i = 1:length(Liste_fich)
    disp([Par Liste_fich{i}]);
end



%% BROUILLON

path = '.' ; % répertoire courant
D = dir(path) ; % récupère un tableau de structure
D = D(~cell2mat({D(:).isdir})) ; % filter pour ne garder que les noms de fichiers
Liste = {D(:).name} ; % transformer en un tableau de cellules texte



%%   AUTRE PROGRAMME 


%----------------------------------
% Fabien Baillon et Jean-Louis Dirion - Nov.2014
% Exemple de traitement par lot
%
clear all
clc
%
experience = 'ManipsLS5';
extension ='csv';
 
%----------------------------------
% Récupération de la liste des fichiers de données présents dans le dossier
% d'expérience
filelist = dir([experience,'/*',extension]);
nfiles = length(filelist);
%
% Pour chaque fichier, on récupère les infos (dates, heures) et les données
% les données sont stockées dans des tableaux RAIES et DATAS
LEGENDS = [];
 
for ifile = 1:nfiles
    disp(['Traitement du fichier n° ',sprintf('%d',ifile)])
    probname = filelist(ifile).name;
    probdate = filelist(ifile).date;
    tmp=csvread([experience,'/',probname]);
%   
    if ~exist('DATAS')
        RAIES = tmp(:,1);
        DATAS = tmp(:,2);
    else
        DATAS = [DATAS tmp(:,2)];
    end
%    
    LEGENDS = [LEGENDS; [experience, ' : ',probdate]];
end
 
 
%----------------------------------
%
plot(RAIES,DATAS);
legend(LEGENDS);
xlabel('\lambda (cm^{-1})');
ylabel('Intensité I (UA)');
title('Spectres RAMAN');


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  ESSAI 


%----------------------------------
% Fabien Baillon et Jean-Louis Dirion - Nov.2014
% Exemple de traitement par lot
%
clear all
clc
%

[F,Par]=uigetfile('*.txt','FICHIER A OUVRIR');
filename=[Par F];
experience = Par;
extension ='txt';
 
%----------------------------------
% Récupération de la liste des fichiers de données présents dans le dossier
% d'expérience
filelist = dir([experience,'/*',extension]);
nfiles = length(filelist);
%
% Pour chaque fichier, on récupère les infos (dates, heures) et les données
% les données sont stockées dans des tableaux RAIES et DATAS
LEGENDS = [];
LEGENDS_n = [];
 
for ifile = 1:nfiles
    %disp(['Traitement du fichier n° ',sprintf('%d',ifile)])
    probname = filelist(ifile).name(1:end-4);
    probdate = filelist(ifile).date;
    %tmp=csvread([experience,'/',probname]);
%   
    
%    
    LEGENDS_n = [LEGENDS_n; [experience, '/',probname]];
    LEGENDS = [LEGENDS; [experience, ' : ',probdate]];
end
 
 
%----------------------------------
%
plot(RAIES,DATAS);
legend(LEGENDS);
xlabel('\lambda (cm^{-1})');
ylabel('Intensité I (UA)');
title('Spectres RAMAN');


for ifile = 1:nfiles
    disp(filelist(ifile).name)
    %disp([num2str(filelist(ifile).name),'%5.2f \t')])
    %probname = filelist(ifile).name;
end






%%
%--------------------------


file_dossier='/Users/sharscoat/Documents/SHARSCOAT/Master_IRIV/Capnographie/0_problemes/209786.txt'
[F,Par]=uigetfile(file_dossier,'FICHIER A OUVRIR');

%[F,Par]=uigetfile('*.txt','FICHIER A OUVRIR');
filename=[Par F];

 
FILE=fopen(filename);
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

