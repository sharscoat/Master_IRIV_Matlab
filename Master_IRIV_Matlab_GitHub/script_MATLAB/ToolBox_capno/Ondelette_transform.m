%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                                                               %%%%
%%%%                       Script de Travail                       %%%%
%%%%                     Transofmé en Ondelette                    %%%%
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



%% sélection du fichier puis ouverture, conversion et récupération des données

%%% Ondelette n°1 - NORMAL
[F,Par]=uigetfile('*.txt','FICHIER A OUVRIR','/Users/sharscoat/Documents/SHARSCOAT/Master_IRIV/Capnographie/wave_form/WaveForm_NORMAL.txt');
filename=[Par F];

FILE=fopen(filename);
Tw1=textscan(FILE,'%s');
fclose(FILE);
Tw1=Tw1{1};
Tw1=strrep(Tw1,',','.');
yw1=cellfun(@str2num,(Tw1(1:end)));
tw1=[0:1:length(yw1)-1];
plot(tw1,yw1)



%%% Ondelette n°2 - BPCO & EMPHYSEMATEUX
[F,Par]=uigetfile('*.txt','FICHIER A OUVRIR','/Users/sharscoat/Documents/SHARSCOAT/Master_IRIV/Capnographie/wave_form/WaveForm_BPCO_EMPHYSEME.txt');
filename=[Par F];

FILE=fopen(filename);
Tw2=textscan(FILE,'%s');
fclose(FILE);
Tw2=Tw2{1};
Tw2=strrep(Tw2,',','.');
yw2=cellfun(@str2num,(Tw2(1:end)));
tw2=[0:1:length(yw2)-1];
plot(tw2,yw2)



%%% Ondelette n°3 - IDENTIFIER
[F,Par]=uigetfile('*.txt','FICHIER A OUVRIR','/Users/sharscoat/Documents/SHARSCOAT/Master_IRIV/Capnographie/wave_form/WaveForm_IDENTIFIER.txt');
filename=[Par F];

FILE=fopen(filename);
Tw3=textscan(FILE,'%s');
fclose(FILE);
Tw3=Tw3{1};
Tw3=strrep(Tw3,',','.');
yw3=cellfun(@str2num,(Tw3(1:end)));
tw3=[0:1:length(yw3)-1];
plot(tw3,yw3)

yw3 = yw3-min(yw3);
inv_yw3 = yw3*(-1);
yw3_p = [yw3(1:132); inv_yw3(28:end)];

tw3_p=[0:1:length(yw3_p)-1];

WaveForm_IDENTIFIER = yw3_p;

filename2=[ '/Users/sharscoat/Documents/SHARSCOAT/Master_IRIV/Capnographie/wave_form/WaveletForm_IDENTIFIER.txt'];
FILE2=fopen(filename2,'w');
fprintf(FILE2,'%4.4f\n',WaveForm_IDENTIFIER);
fclose(FILE2);



%%% Signal n°1.
[F,Par]=uigetfile('*.txt','FICHIER A OUVRIR','/Users/sharscoat/Documents/SHARSCOAT/Master_IRIV/Capnographie/wave_form/WaveForm_NORMAL_25.txt');
filename=[Par F];

FILE=fopen(filename);
S1=textscan(FILE,'%s');
fclose(FILE);
S1=S1{1};
S1=strrep(S1,',','.');
sy1=cellfun(@str2num,(S1(1:end)));
st1=[0:1:length(sy1)-1];
plot(st1,sy1);








%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Transformé Discrète en Ondelettes    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%load noissin;
%whos
%sig = noissin;
sig = sy1;
longing = length(sig);
figure; plot(sig); print -depsc2 ex2-1
[a,d] = dwt(sig,?db2?);
figure; plot(a); print -depsc2 ex2-2
figure; plot(d); print -depsc2 ex2-3
a1 = upcoef(?a?,a,?db2?,1,longsig);
d1 = upcoef(?d?,d,?db2?,1,longsig);
figure; plot(a1); print -depsc2 ex2-4
figure; plot(d1); print -depsc2 ex2-5
a0 = idwt(a,d,?db2?,longsig);
figure; plot(a0); print -depsc2 ex2-6
erreur = sig-a0;
figure; plot(erreur); print -depsc2 ex2-7



L=12;
W=FWT_PO(sig,L,qmf);

qmf=MakeONFilter(?Haar?);
qmf2=MakeONFilter(?Daubechies?,4);






%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Intercorrélation : function "xcorr()"  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% signal de base
%s1 = sy1;   % signal  
%s2 = yw1;   % motif de base

%s1 = yw1;
%s2 = yw1;

s1 = yw1;
s2 = yw2;

%s1 = yw1;
%s2 = yw3;


% temps et fréquence
nue = 50;                   % fréquence d'échantillonnage
te = 1/nue;                 % période d'échantillonnage
nacq = length(sy1);           % nombre d'échantillons
tacq = (nacq-1)*te;         % durée d'acquisition
t = [0:te:tacq];            % vecteur des temps d'acquisition
dnu = 1/tacq;               % résolution fréquentielle
numax = nue/2;              % fréquence maximale
nus = [-numax:dnu:numax];   % vecteur des fréquences

Fs = nue;   % fréquence d'aquisition

t1 = (0:length(s1)-1)/Fs;
t2 = (0:length(s2)-1)/Fs;

% graphique
subplot(2,1,1)
plot(t1,s1)
title('s_1')

subplot(2,1,2)
plot(t2,s2)
title('s_2')
xlabel('Time (s)')

% Intercorrélation 
[acor,lag] = xcorr(s2,s1);

[~,I] = max(abs(acor));
lagDiff = lag(I)


timeDiff = lagDiff/Fs


figure
plot(lag,acor)
a3 = gca;
a3.XTick = sort([-3000:1000:3000 lagDiff]);


