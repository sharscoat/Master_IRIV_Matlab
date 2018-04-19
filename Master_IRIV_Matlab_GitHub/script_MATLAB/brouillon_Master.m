%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                                                               %%%%
%%%%                       Script de Travail                       %%%%
%%%%                          CAPNOGRAPHY                          %%%%
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



%% EN fait, on va essayer de trouver les différents cycles sans findpeaks

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

%% indexation terminée
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Extration des différents cycles : CAPNOGRAMMES

%hwb=waitbar(0,'identification des périodes');
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
%close(hwb);

p=input('quel capnogramme : ');
plot(R(p).tp,R(p).yp,'g')

% les différents cycles ou Capnogrammes sont mis dans un "Objet" ou "Liste"
% ou "Dictionnaire" ? appelé : "R"
% avec pour composante "tp" et "yp"
% avec un nombre de sous-calsse : "p"  (indexé de 1 à p)

hold on
for p=1:NP
    plot(R(p).yp)
end



%% Capnographe Pathologique typique

%%%%%%%%%%%%%%%%%%%%%  BPCO - EMPHYSEME %%%%%%%%%%%%%%%%%%%%%
%%%% Wave Form de patient pathologique typique BPCO - EMPHYSEME
capno_patho = y(B(7):B(8));
temps_patho = t(B(7):B(8));
% ou 
t3 = R(7).tp;
y3 = R(7).yp;

plot(temps_patho,capno_patho);
% les capnogrammes 3 et 7 sont de bonne qualité pour une forme de courbe
% pahtologique

wave_form_5 = [capno_patho; capno_patho; capno_patho; capno_patho; capno_patho];
wave_form_25 = [wave_form_5; wave_form_5; wave_form_5; wave_form_5; wave_form_5 ];
temps_form_5 = [0:1:length(wave_form_5)-1];

plot(temps_form_5,wave_form_5);

filename2=[ '/Users/sharscoat/Documents/SHARSCOAT/Master_IRIV/Capnographie/wave_form/WaveForm_BPCO_EMPHYSEME.txt'];
FILE2=fopen(filename2,'w');
fprintf(FILE2,'%4.4f\n',y3);
fclose(FILE2);


filename2=[ '/Users/sharscoat/Documents/SHARSCOAT/Master_IRIV/Capnographie/wave_form/WaveForm_BPCO_EMPHYSEME_25.txt'];
FILE2=fopen(filename2,'w');
fprintf(FILE2,'%4.4f\n',wave_form_25);
fclose(FILE2);


%%%%%%%%%%%%%%%%%%%%%  NORMAL %%%%%%%%%%%%%%%%%%%%%
%%%% Wave Form de patient normal
capno_patho = y(B(2):B(3));
temps_patho = t(B(2):B(3));
% ou 
t2 = R(2).tp;
y2 = R(2).yp;
plot(t2,y2);


wave_form_5 = [y2; y2; y2; y2; y2];
wave_form_25 = [wave_form_5; wave_form_5; wave_form_5; wave_form_5; wave_form_5];
temps_form_5 = [0:1:length(wave_form_5)-1];

plot(temps_form_5,wave_form_5);


filename2=[ '/Users/sharscoat/Documents/SHARSCOAT/Master_IRIV/Capnographie/wave_form/WaveForm_NORMAL.txt'];
FILE2=fopen(filename2,'w');
fprintf(FILE2,'%4.4f\n',y2);
fclose(FILE2);


filename2=[ '/Users/sharscoat/Documents/SHARSCOAT/Master_IRIV/Capnographie/wave_form/WaveForm_NORMAL_25.txt'];
FILE2=fopen(filename2,'w');
fprintf(FILE2,'%4.4f\n',wave_form_25);
fclose(FILE2);





%%%%%%%%%%%%%%%%%%%%%  BPCO %%%%%%%%%%%%%%%%%%%%%
%%%% Wave Form de patient pathologique type BPCO
capno_patho = y(B(4):B(5));
temps_patho = t(B(4):B(5));
% ou 
t4 = R(4).tp;
y4 = R(4).yp;
plot(t4,y4);


wave_form_5 = [y4; y4; y4; y4; y4];
wave_form_25 = [wave_form_5; wave_form_5; wave_form_5; wave_form_5; wave_form_5];
temps_form_5 = [0:1:length(wave_form_5)-1];

plot(temps_form_5,wave_form_5);


filename2=[ '/Users/sharscoat/Documents/SHARSCOAT/Master_IRIV/Capnographie/wave_form/WaveForm_BPCO.txt'];
FILE2=fopen(filename2,'w');
fprintf(FILE2,'%4.4f\n',y4);
fclose(FILE2);


filename2=[ '/Users/sharscoat/Documents/SHARSCOAT/Master_IRIV/Capnographie/wave_form/WaveForm_BPCO_25.txt'];
FILE2=fopen(filename2,'w');
fprintf(FILE2,'%4.4f\n',wave_form_25);
fclose(FILE2);





%%%%%%%%%%%%%%%%%%%%%  A IDENTIFIER %%%%%%%%%%%%%%%%%%%%%
%%%% Wave Form de patient A Identifier
capno_patho = y(B(8):B(9));
temps_patho = t(B(8):B(9));
% ou 
t8 = R(8).tp;
y8 = R(8).yp;
plot(t8,y8);


wave_form_5 = [y8; y8; y8; y8; y8];
wave_form_25 = [wave_form_5; wave_form_5; wave_form_5; wave_form_5; wave_form_5];
temps_form_5 = [0:1:length(wave_form_5)-1];

plot(temps_form_5,wave_form_5);


filename2=[ '/Users/sharscoat/Documents/SHARSCOAT/Master_IRIV/Capnographie/wave_form/WaveForm_IDENTIFIER.txt'];
FILE2=fopen(filename2,'w');
fprintf(FILE2,'%4.4f\n',y8);
fclose(FILE2);


filename2=[ '/Users/sharscoat/Documents/SHARSCOAT/Master_IRIV/Capnographie/wave_form/WaveForm_IDENTIFIER_25.txt'];
FILE2=fopen(filename2,'w');
fprintf(FILE2,'%4.4f\n',wave_form_25);
fclose(FILE2);





%%
%--------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  Lecture d'un Fichier PDF %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



filename_pdf=[ '/Users/sharscoat/Documents/SHARSCOAT/Master_IRIV/Capnographie/wave_form/essai_EFR.pdf'];
FILE_pdf=fopen(filename_pdf,'r');
contenu_pdf=fscanf(FILE_pdf,'%s');


%%%%%%%%%%%%%


filename=[ '/Users/sharscoat/Documents/SHARSCOAT/Master_IRIV/Capnographie/wave_form/essai_EFR.txt'];
FILE=fopen(filename,'r');
T=textscan(FILE,'%s');
fclose(FILE);
T=T{1};
T{60:70}







%%
%--------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%


%% sélection du fichier

[F,Par]=uigetfile('*.txt','FICHIER A OUVRIR');
filename=[Par F];

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

y = y - min(y);
y = y(19:end);


% inversion du vecteur y
y_inv= y(length(y):-1:1);


% inversion de la courbe sur l'axe des y
y_inv = -1*y_inv;



% concaténation des deux vecteurs
yy = [y_inv; y];

plot(yy);ylim([-6 6]); grid on;




