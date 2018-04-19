%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                                                               %%%%
%%%%                       Script de Travail                       %%%%
%%%%                 "Visualisation de Capnographie"               %%%%
%%%%                       "FILTRE PASSE BAS"                      %%%%
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



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%   FILTRE PASSE BAS  %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



nue = 50;     % fréquence d'échantillonnage
te = 1/nue;      % période d'échantillonnage
nacq = length(y);	% nombre d'échantillons
tacq = (nacq-1)*te;	% durée d'acquisition
t = [0:te:tacq];       % vecteur des temps d'acquisition
dnu = 1/tacq;     % résolution fréquentielle
numax = nue/2;   % fréquence maximale
nus = [-numax:dnu:numax];     % vecteur des fréquences


% calcul du spectre
sy = fftshift(fft(y));

% partie réel
Ry = real(sy);
% partie imaginaire
Iy = imag(sy);
% module 
My = abs(sy);


% représentation graphique
figure(1);
subplot(2,2,1)
plot(t(1:500),y(1:500)); % signal
    title('Signal');
    set(get(gca,'XLabel'),'String','Temps (s)','FontName','Arial','FontSize',9)
    set(get(gca,'YLabel'),'String','Amplitude (V)','FontName','Arial','FontSize',9)
subplot(2,2,2)
plot(nus,Ry);	% partie réelle du spectre
    title('partie réelle');
    set(get(gca,'XLabel'),'String','Temps (s)','FontName','Arial','FontSize',9)
    set(get(gca,'YLabel'),'String','Fréquence (Hz)','FontName','Arial','FontSize',9)
subplot(2,2,3)
plot(nus,Iy);	% partie imaginaire du spectre
    title('partie imaginaire');
    set(get(gca,'XLabel'),'String','Temps (s)','FontName','Arial','FontSize',9)
    set(get(gca,'YLabel'),'String','Fréquence (Hz)','FontName','Arial','FontSize',9)
subplot(2,2,4)
plot(nus,My);	% module du spectre
    title('partie module');
    set(get(gca,'XLabel'),'String','Temps (s)','FontName','Arial','FontSize',9)
    set(get(gca,'YLabel'),'String','Fréquence (Hz)','FontName','Arial','FontSize',9)




%%%%%%%%%%%%%%%%%%%%%%
%%%% filtre passe bas
%%%%%%%%%%%%%%%%%%%%%

sy = fftshift(fft(y));
sy_filt_pb = sy;
ft_pb = [1:1:length(y)];
centre_x = floor(length(y)/2);     
N = 100;
sy_fpb = sy_filt_pb;
sy_fpb(1:centre_x - N) = 0;
sy_fpb(centre_x + N : end) = 0;
% transformation inverse
y_fpb = ifft(ifftshift(sy_fpb));
% module 
Myf = abs(sy_fpb);
plot(nus,Myf);


%%%%%  affichage des 2 courbes : native & filtre PB  %%%%%
Par=get(0,'ScreenSize');
h=figure;
set(h,'Position',Par);
set(h,'PaperOrientation','landscape');
% set(h,'Position',[50 50 1200 800]);
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
subplot(2,1,1)
plot(t,y,'b', 'MarkerFaceColor','b')
grid
hold on
subplot(2,1,2)
plot(t,y_fpb,'r', 'MarkerFaceColor','r')
grid







%%%%%%%%%%%%%%%%%%%%%%
%%%% filtre moyenneur
%%%%%%%%%%%%%%%%%%%%%

sy = fftshift(fft(y));
sy_filt_pb = sy;
ft_pb = [1:1:length(y)];
centre_x = floor(length(y)/2);     
N = 100;
sy_fpb = sy_filt_pb;
sy_fpb(1:centre_x - N) = 0;
sy_fpb(centre_x + N : end) = 0;
% transformation inverse
y_fpb = ifft(ifftshift(sy_fpb));
% module 
Myf = abs(sy_fpb);
plot(nus,Myf);


%%%%%  affichage des 2 courbes : native & filtre PB  %%%%%
Par=get(0,'ScreenSize');
h=figure;
set(h,'Position',Par);
set(h,'PaperOrientation','landscape');
% set(h,'Position',[50 50 1200 800]);
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
subplot(2,1,1)
plot(t,y,'b', 'MarkerFaceColor','b')
grid
hold on
subplot(2,1,2)
plot(t,y_fpb,'r', 'MarkerFaceColor','r')
grid


















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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











%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%   FILTRE PASSE BAS  %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



nue = 50;     % fréquence d'échantillonnage
te = 1/nue;      % période d'échantillonnage
nacq = length(y);	% nombre d'échantillons
tacq = (nacq-1)*te;	% durée d'acquisition
t = [0:te:tacq];       % vecteur des temps d'acquisition
dnu = 1/tacq;     % résolution fréquentielle
numax = nue/2;   % fréquence maximale
nus = [-numax:dnu:numax];     % vecteur des fréquences


% calcul du spectre
sy = fftshift(fft(y));

% partie réel
Ry = real(sy);
% partie imaginaire
Iy = imag(sy);
% module 
My = abs(sy);


% représentation graphique
figure(1);
subplot(2,2,1)
plot(t(1:500),y(1:500)); % signal
    title('Signal');
    set(get(gca,'XLabel'),'String','Temps (s)','FontName','Arial','FontSize',9)
    set(get(gca,'YLabel'),'String','Amplitude (V)','FontName','Arial','FontSize',9)
subplot(2,2,2)
plot(nus,Ry);	% partie réelle du spectre
    title('partie réelle');
    set(get(gca,'XLabel'),'String','Temps (s)','FontName','Arial','FontSize',9)
    set(get(gca,'YLabel'),'String','Fréquence (Hz)','FontName','Arial','FontSize',9)
subplot(2,2,3)
plot(nus,Iy);	% partie imaginaire du spectre
    title('partie imaginaire');
    set(get(gca,'XLabel'),'String','Temps (s)','FontName','Arial','FontSize',9)
    set(get(gca,'YLabel'),'String','Fréquence (Hz)','FontName','Arial','FontSize',9)
subplot(2,2,4)
plot(nus,My);	% module du spectre
    title('partie module');
    set(get(gca,'XLabel'),'String','Temps (s)','FontName','Arial','FontSize',9)
    set(get(gca,'YLabel'),'String','Fréquence (Hz)','FontName','Arial','FontSize',9)





%%%% filtre passe bas

sy = fftshift(fft(y));
sy_filt_pb = sy;
ft_pb = [1:1:length(y)];
centre_x = floor(length(y)/2);     
N = 50;
sy_fpb = sy_filt_pb;
sy_fpb(1:centre_x - N) = 0;
sy_fpb(centre_x + N : end) = 0;
% transformation inverse
y_fpb = ifft(ifftshift(sy_fpb));

plot(t,y_fpb);

%%%%%  affichage des 2 courbes : native & filtre PB  %%%%%
Par=get(0,'ScreenSize');
h=figure;
set(h,'Position',Par);
set(h,'PaperOrientation','landscape');
% set(h,'Position',[50 50 1200 800]);
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
subplot(2,1,1)
plot(t,y,'b', 'MarkerFaceColor','b')
grid
hold on
subplot(2,1,2)
plot(t,y_fpb,'r', 'MarkerFaceColor','r')
grid







% partie réel
Ryf = real(sy_fpb);
% partie imaginaire
Iyf = imag(sy_fpb);
% module 
Myf = abs(sy_fpb);

% plot(nus,Myf);

















