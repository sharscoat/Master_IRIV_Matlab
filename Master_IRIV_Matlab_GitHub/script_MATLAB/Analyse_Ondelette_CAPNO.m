%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                                                               %%%%
%%%%                  Script d'Analyse en Ondelette                %%%%
%%%%                          CAPNOGRAPHY                          %%%%
%%%%                                                               %%%%
%%%%                       S�bastien Harscoat                      %%%%
%%%%                                                               %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% initialisation
clc
clear
close all



%% s�lection du fichier

[F,Par]=uigetfile('*.txt','FICHIER A OUVRIR');
filename=[Par F];
% filename2=[ filename(1:end-4) '_.csv'];
% filename3=[ filename(1:end-4) '.pdf'];
filename4=[ filename(1:end-4) '_Res.csv'];
% '/Users/sharscoat/Documents/SHARSCOAT/Master_IRIV/Capnographie/MatLab/scripts/10581GUHMANN_EMPHYSEME.txt'

% NOM du FICHIER
n_fich=filename(end-15:end);
debut=strfind(n_fich,'/')+1;
fin=strfind(n_fich,'.')-1;
nom_fich=n_fich(debut:fin);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% initialisation
clc
clear
close all


%% r�cuperer les noms des fichiers dans un dossier

%[F,Par]=uigetfile('*.txt','FICHIER A OUVRIR');
%filename=[Par F];

file_dossier='/Users/sharscoat/Documents/SHARSCOAT/Master_IRIV/Capnographie/essaie/10581.txt'
[F,Par]=uigetfile(file_dossier,'FICHIER A OUVRIR');

path = Par ; % r�pertoire courant
D = dir(path); % r�cup�re un tableau de structure
D = D(~cell2mat({D(:).isdir})); % filter pour ne garder que les noms de fichiers
Liste = {D(:).name}; % transformer en un tableau de cellules texte
Liste = Liste';
Liste_fich = Liste(2:length(Liste));

% tous les chemin pour chaque fichier :

for i = 1:length(Liste_fich)
    liste_path{i}=[Par Liste_fich{i}];
end



%% Analyse de tous les fichiers

for j = 1:length(liste_path)
    %% ouverture, conversion et r�cup�ration des donn�es
    %filename = liste_path{j};
    dir(liste_path{j})
    FILE=fopen(liste_path{j});
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

    %%% la  portion utile de la courbe est d�sormais choisie



    %% EN fait, on va essayer de trouver les diff�rents cycles sans findpeaks

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

    %% indexation termin�e
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %% Extration des diff�rents cycles : CAPNOGRAMMES

    %hwb=waitbar(0,'identification des p�riodes');
    for p=1:NP
        %waitbar(p/NP,hwb);
        indp=B(p):B(p+1);
        tp=t(indp);
        yp=y(indp);
        R(p).tp=tp;
        R(p).yp=yp;

        MM=max(yp);
        mm=min(yp);

    end

    % plot(R(p).tp,R(p).yp,'g')



    %% constitution de la base de donn�es des coefficients des SCALOGRAM

    for p = 1:NP
            % optimisation de la courbe y
                % �talonage � 0
                % normalisation
                % 
            yy = R(p).yp;
            yy = yy - min(yy);
            yy = yy/max(yy);
            ZZ = find(yy > 0.1,1,'first');
            yy = yy(ZZ-10:end);
            %y = R(p).yp;

            % SCALOGRAM
            data=yy;
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
            %resultat=[fliplr(scalogramme)]';
            % normamisation du scalogramme vis � vis de l'�chantillonage 
            resultat_1=[fliplr(scalogramme)]';
            [resultat] = resultat_1/sqrt(N0);


        An2(p) = max(resultat(2,:));
        An3(p) = max(resultat(3,:));
        An4(p) = max(resultat(4,:));
        An5(p) = max(resultat(5,:));

        An(1,p) = max(resultat(1,:));
        An(2,p) = max(resultat(2,:));
        An(3,p) = max(resultat(3,:));
        An(4,p) = max(resultat(4,:));
        An(5,p) = max(resultat(5,:));
        An(6,p) = max(resultat(6,:));
        [An] = An;
    end
    
    
    %base{j,1} = An;
    base{j,1} = An(1,:);
    base{j,2} = An(2,:);
    base{j,3} = An(3,:);
    base{j,4} = An(4,:);
    base{j,5} = An(5,:);
    base{j,6} = An(6,:);
    base{j,7} = Liste_fich{j}(1:end-4);
    
end


%% Ectiture dans un fichier ".csv"

%[F,Par]=uigetfile('*.csv','FICHIER A OUVRIR');
%filename_res=[Par F];

filename_res='/Users/sharscoat/Documents/SHARSCOAT/R_project/Master_IRIV/resultats_scalogram.csv'


% nettoyage du fichier avant �criture
FILE = fopen(filename_res,'w');
fprintf(FILE,'\n');
fclose(FILE);

% Ecriture dans les fichier "resultats_scalogram.csv"
FILE = fopen(filename_res,'w');
for j = 1:length(liste_path)
    fprintf(FILE,'\n');
    fprintf(FILE,base{j,7});
    fprintf(FILE,'\n');
    fprintf(FILE,'%f;',base{j,1});
    fprintf(FILE,'\n');
    fprintf(FILE,'%f;',base{j,2});
    fprintf(FILE,'\n');
    fprintf(FILE,'%f;',base{j,3});
    fprintf(FILE,'\n');
    fprintf(FILE,'%f;',base{j,4});
    fprintf(FILE,'\n');
    fprintf(FILE,'%f;',base{j,5});
    fprintf(FILE,'\n');
    fprintf(FILE,'%f;',base{j,6});
end
fclose(FILE);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%




base{j,1} = An(1,:);
base{j,2} = An(2,:);
base{j,3} = An(3,:);
base{j,4} = An(4,:);
base{j,5} = An(5,:);
base{j,6} = An(6,:);
base{j,7} = Liste_fich{j}(1:end-4);



figure;
hold on;
plot(An2,'m.-'); plot(An3,'g.-'); ;plot(An4,'r.-'); plot(An5,'b.-');

figure;
plot(An');


%% Ectiture dans un fichier ".csv"


FILE = fopen(filename4,'w');
%fprintf(FILE,'\n');
fprintf(FILE,nom_fich);
fprintf(FILE,'\n');
fprintf(FILE,'%f;',An(1,:));
fprintf(FILE,'\n');
fprintf(FILE,'%f;',An(2,:));
fprintf(FILE,'\n');
fprintf(FILE,'%f;',An(3,:));
fprintf(FILE,'\n');
fprintf(FILE,'%f;',An(4,:));
fprintf(FILE,'\n');
fprintf(FILE,'%f;',An(5,:));
fprintf(FILE,'\n');
fprintf(FILE,'%f;',An(6,:));
fclose(FILE);






%%  BROUILLON
%-----------------------

% ouverture du fichier lecture et r��criture

FILE = fopen(filename4,'w');
%fprintf(FILE,'%f;;',An);
fprintf(FILE,'\n');
fprintf(FILE,nom_fich);
fprintf(FILE,'\n');
fprintf(FILE,'%f;',An2);
fprintf(FILE,'\n');
fprintf(FILE,'%f;',An3);
fprintf(FILE,'\n');
fprintf(FILE,'%f;',An4);
fprintf(FILE,'\n');
fprintf(FILE,'%f;',An5);
fclose(FILE);


FILE=fopen(filename4);
Tn=textscan(FILE,'%s');
fclose(FILE);
Tn=Tn{1};
T=strrep(T,',','.');
T(1:10)=[];
T(end-2:end)=[];
c=T(2:2:end);
y=cellfun(@str2num,(c(1:end)));





