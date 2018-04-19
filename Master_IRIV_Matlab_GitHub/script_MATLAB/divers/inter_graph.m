%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INTERFACE GRAPHIQUE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DEBUT DE LA FONCTION PRINCIPALE%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function inter_graph

% Définition de nCompteur et handles comme variables globales dans chaque fonction et sous-fonction
% nCompteur : valeur courante du compteur (scalaire)
% handles : identifiants des objets graphiques (vecteur)
global nCompteur handles

% Initialisation de la variable représentant la valeur courante du compteur nCompteur à 0
nCompteur=0;

% Création de l'objet Figure
handles(1)=figure('units','pixels',...
    'position',[250 250 500 500],...
    'color',[0.925 0.913 0.687],...
    'numbertitle','off',...
    'name','[GUI] Utilisation des variables globales',...
    'menubar','none',...
    'tag','interface');

% Création de l'objet Uicontrol Pushbutton -
handles(2)=uicontrol('style','pushbutton',...
    'units','normalized',...
    'position',[0.1 0.1 0.1 0.05],...
    'string','-',...    
    'callback',@retrancher,...
    'tag','bouton-');

% Création de l'objet Uicontrol Pushbutton +
handles(3)=uicontrol('style','pushbutton',...
    'units','normalized',...
    'position',[0.3 0.1 0.1 0.05],...
    'string','+',...    
    'callback',@ajouter,...
    'tag','bouton+');

% Création de l'objet Uicontrol Text résultat
handles(4)=uicontrol('style','text',...
    'units','normalized',...
    'position',[0.1 0.2 0.3 0.05],...
    'string','0',...
    'tag','resultat');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FIN DE LA FONCTION PRINCIPALE%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DEBUT DE LA SOUS-FONCTION RETRANCHER%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function retrancher(obj,event)

% Définition de nCompteur et handles comme variables globales dans chaque fonction et sous-fonction
% nCompteur : valeur courante du compteur (scalaire)
% handles : identifiants des objets graphiques (vecteur)
global nCompteur handles

% Diminution de la valeur de nCompteur
nCompteur=nCompteur-1;
 
% Actualisation de la propriété String de l'objet Uicontrol Text résultat
set(handles(4),'string',num2str(nCompteur));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FIN DE LA SOUS-FONCTION RETRANCHER%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DEBUT DE LA SOUS-FONCTION AJOUTER%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ajouter(obj,event)

% Définition de nCompteur et handles comme variables globales dans chaque fonction et sous-fonction
% nCompteur : valeur courante du compteur (scalaire)
% handles : identifiants des objets graphiques (vecteur)
global nCompteur handles

% Augmentation de la valeur de nCompteur
nCompteur=nCompteur+1;

% Actualisation de la propriété String de l'objet Uicontrol Text résultat
set(handles(4),'string',num2str(nCompteur));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FIN DE LA SOUS-FONCTION AJOUTER%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








%%%%%%%%%%%%%%%%%

%handles = guidata(handles_de_la_fenetre);
 
%number = str2num(get(handles.edit, 'String'));