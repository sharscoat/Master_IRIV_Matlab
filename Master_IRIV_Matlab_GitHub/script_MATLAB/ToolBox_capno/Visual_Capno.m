%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                                                               %%%%
%%%%                       Script de Travail                       %%%%
%%%%                 "Visualisation de Capnographie"               %%%%
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

[F,Par]=uigetfile('*.txt','FICHIER A OUVRIR','/Users/sharscoat/Documents/SHARSCOAT/Master_IRIV/Capnographie/CAPNO_ref_propre/');
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

%p=input('quel capnogramme : ');
%plot(R(p).tp,R(p).yp,'g')







%%
%%%%%%%%%%%%%%  AFFICHAGE %%%%%%%%%%%%%%

%% affichage
Par=get(0,'ScreenSize');
h=figure;
set(h,'Position',Par);
set(h,'PaperOrientation','landscape');
% set(h,'Position',[50 50 1200 800]);
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);

%p=2;

%% cycle d'affichage
for p=1:NP
indp=B(p):B(p+1);
tp=t(indp);
yp=y(indp);


subplot(3,1,1)
plot(t,y,'r', 'MarkerFaceColor','r')
grid
hold on
plot(tp,yp,'o','Color',[.5 0 0])
plot(t(B),y(B),'sk','MarkerFaceColor','k', 'MarkerSize',8)
%hold on

subplot(3,2,[3 5])
plot(tp,yp,'r')
grid
%hold on

subplot(3,2,6)
cla
axis off

% pause
waitfor(msgbox('NEXT','APPUYER POUR CONTINUER'))
end











