% ajustement SH

%% initialisation
clc
clear
%% 
close all


DelaiPente=0;
while DelaiPente<1
    txt=cell2mat(inputdlg('Nombre de points pour l''ajustement linéaire','Points ?'));
    answer = ceil(str2num(txt));
    if ~isempty(answer) && answer>0 && answer<100
        DelaiPente=answer;
    end
end

Tfonc='c+((d-c)-a*(t1-t)).*(t>t0).*(t<=t1).*(1-exp(-(t-t0)/T1))+ (d-c)*exp(-(t-t1)/T2).*(t>t1)';

%% sélection du fichier

[F,Par]=uigetfile('*.txt','FICHIER CSV A OUVRIR');
filename=[Par F];
filename2=[ filename(1:end-4) '_.csv'];
filename3=[ filename(1:end-4) '.pdf'];
filename4=[ filename(1:end-4) '_Res.csv'];

%% ouverture, conversion et récupération des données
FILE=fopen(filename);
T=textscan(FILE,'%s');
fclose(FILE);
T=T{1};
T=strrep(T,',','.');
K={sprintf('%s\n', T{:})};
L=cell2mat(K);

FILE2=fopen(filename2,'w');
fprintf(FILE2,'%s',L);
fclose(FILE2);

FILE2=fopen(filename2);
U=textscan(FILE2,'%s',2);
V=textscan(FILE2,'%s%f','Delimiter',',');
fclose(FILE2);
clear U
a=V{1};
clear V
t=cell2mat(a(1:2:end));
m=t(:,1:2);
t(:,1:3)=[];
t(:,end)=[];
t=str2num(t)+60*str2num(m); %#ok<*ST2NM>
y=cellfun(@str2num,(a(2:2:end)));

% % % le problème ici est le signe moins éventuel 
% % % il faudrait ajouter un espace avant s'il n'y a pas de signe moins....
% % %     ou changer la fonction
% % % y=str2num(cell2mat(a(2:2:end)));
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

%% la  portion utile de la courbe est désormais choisie
%% il va falloir désormais sélectionner le cycle choisi manuellement
%% ou bien trouver tous les cycles (pour ceci, il faut la fonction findpeaks du Toolbox Signal Processing

%% pour l'instant, on prend le premier pic


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
%% calcul des optimisations sur les différents cycles
hwb=waitbar(0,'identification des périodes');
for p=1:NP
    waitbar(p/NP,hwb);
    indp=B(p):B(p+1);
    tp=t(indp);
    yp=y(indp);
    R(p).tp=tp;
    R(p).yp=yp;
    
    MM=max(yp);
    mm=min(yp);
    
    %% dérivées
    t1=diff(tp);
    y1=diff(yp)./t1;
    y2=diff(y1)./t1(2:end);

    %% trouver t0
    Ma=max(y1);
    t0_=find(y1>Ma/10, 1, 'first');    % index
    t0=tp(t0_);                        % valeur (en s)

    %% touver t1
    Mi=min(y1);
    t1_=find(y1<Mi/10, 1,'first');
    t1=tp(t1_);

    %% trouver c
    c=mean(yp(1:t0_));

    %% trouver d
    d=MM;

    %% trouver a
    tt0_=floor(t1_-(t1_-t0_)/3);
    tt1_=t1_-3;
    a=median(diff(yp(tt0_:tt1_))./diff(tp(tt0_:tt1_)));
%     std(diff(yp(tt0_:tt1_))./diff(tp(tt0_:tt1_)))    
    
    %% estimer T1
    T1=(t1-t0)/2;

    %% estimer T2
    T2=(max(tp)-t1)/2;

    %% AJUSTEMENT
    [xData, yData] = prepareCurveData( tp, yp );
    
    FT = fittype( Tfonc, 'independent', 't' );

    TM=max(tp);

    opts = fitoptions( FT );
    opts.MaxIter = 500;
    coeffnames(FT);

    opts.Display = 'Off';

    opts.StartPoint =   [	T1       T2     a       c       d       t0      t1     ];
    opts.Lower =        [   eps     eps     eps   mm     mm      0       0       ];
    opts.Upper =        [   10*TM   10*TM   9999    MM     MM       TM      TM     ];

    [cfun,gof, output] = fit( xData, yData, FT, opts );
    ci=confint(cfun);
    parms=coeffvalues(cfun);
    T1=parms(1);
    T2=parms(2);
    a =parms(3);
    c =parms(4);
    d =parms(5);
    t0=parms(6);
    t1=parms(7);
    R(p).T1=T1; %#ok<*SAGROW>
    R(p).T2=T2;
    R(p).a =a;
    R(p).c =c;
    R(p).d =d;
    R(p).t0=t0;
    R(p).t1=t1;

    R(p).T1_ic=ci(:,1);
    R(p).T2_ic=ci(:,2);
    a_ic=ci(:,3);
    R(p).a_ic=a_ic;
    R(p).c_ic=ci(:,4);
    R(p).d_ic=ci(:,5);
    R(p).t0_ic=ci(:,6);
    R(p).t1_ic=ci(:,7);


    y_=feval(cfun,tp);
    R(p).y_=y_;
    drt=[min(tp) t1];
    R(p).drt=drt;
    R(p).dr =c + (d-c-a*(t1-drt));
    R(p).dri=c + (d-c-a_ic(1)*(t1-drt));
    R(p).drs=c + (d-c-a_ic(2)*(t1-drt));

    R(p).r_=yp-y_;
    
    R(p).r2=gof.rsquare;
%     h=figure;
%     plot(tp,yp,'or')
%     hold on
%     plot(tp,y_,'k')
%     pause
%     close(h)
end
close(hwb);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% on passe à l'affichage
% choix du cycle médian pour l'affichage initial

%% affichage
Par=get(0,'ScreenSize');
h=figure;
set(h,'Position',Par);
set(h,'PaperOrientation','landscape');
% set(h,'Position',[50 50 1200 800]);
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);


p=2;
%% cycle d'affichage
for p=1:NP
indp=B(p):B(p+1);
tp=t(indp);
yp=y(indp);


subplot(3,1,1)
plot(t,y,'or', 'MarkerFaceColor','r')
grid
hold on
plot(tp,yp,'o','Color',[.5 0 0])
plot(t(B),y(B),'sk','MarkerFaceColor','k', 'MarkerSize',8)


subplot(3,2,[3 5])
plot(tp,yp,'or')
grid
hold on
drt=R(p).drt;
dri=R(p).dri;
drs=R(p).drs;
dr=R(p).dr;


plot(drt,dr,'-','Color',[.3 .3 .5], 'LineWidth', 3) 
plot(drt,dri,'Color',[.5 .5 .5], 'LineWidth', 2) 
plot(drt,drs,'Color',[.5 .5 .7], 'LineWidth', 2) 

t0=R(p).t0;
t1=R(p).t1;
T1=R(p).T1;
T2=R(p).T2;
a=R(p).a;
c=R(p).c;
d=R(p).d;

r2=R(p).r2;

t0_ic(1)=R(p).t0_ic(1);
t1_ic(1)=R(p).t1_ic(1);
T1_ic(1)=R(p).T1_ic(1);
T2_ic(1)=R(p).T2_ic(1);
a_ic(1)=R(p).a_ic(1);
c_ic(1)=R(p).c_ic(1);
d_ic(1)=R(p).d_ic(1);

t0_ic(2)=R(p).t0_ic(2);
t1_ic(2)=R(p).t1_ic(2);
T1_ic(2)=R(p).T1_ic(2);
T2_ic(2)=R(p).T2_ic(2);
a_ic(2)=R(p).a_ic(2);
c_ic(2)=R(p).c_ic(2);
d_ic(2)=R(p).d_ic(2);


y_=R(p).y_;

plot(tp,y_,'k','LineWidth',3)
xlim([min(tp) max(tp)])

subplot(3,2,4)
r_=R(p).r_;
plot(tp,r_,'rx')
xlim([min(tp) max(tp)])

subplot(3,2,6)
cla
axis off

text(0, 1.0, ['Pente (a) = ' num2str(a,'%5.2f \t(')  num2str(a_ic(1),'%5.2f ; ')  num2str(a_ic(2),'%5.2f)')], 'FontSize', 16, 'FontWeight', 'bold') 
text(0, 0.85, ['Offset (c) = ' num2str(c,'%5.2f \t(')  num2str(c_ic(1),'%5.2f ; ')  num2str(c_ic(2),'%5.2f)')], 'FontSize', 14) 
text(0, 0.70, ['Amplitude (d) = ' num2str(d,'%5.2f \t(')  num2str(d_ic(1),'%5.2f ; ')  num2str(d_ic(2),'%5.2f)')], 'FontSize', 14) 
text(0, 0.55, ['Temps montée (t0) = ' num2str(t0,'%5.2f \t(')  num2str(t0_ic(1),'%5.2f ; ')  num2str(t0_ic(2),'%5.2f)')], 'FontSize', 14) 
text(0, 0.40, ['Temps descente (t1) = ' num2str(t1,'%5.2f \t(')  num2str(t1_ic(1),'%5.2f ; ')  num2str(t1_ic(2),'%5.2f)')], 'FontSize', 14) 
text(.7, 0.55, ['délai ( t1 - t0 ) = ' num2str(t1-t0,'%5.2f \t') ], 'FontSize', 14) 
text(0., 0.25, ['Constante temps montée (T1) = ' num2str(T1,'%5.2f \t(')  num2str(T1_ic(1),'%5.2f ; ')  num2str(T1_ic(2),'%5.2f)')], 'FontSize', 14) 
text(0., 0.10, ['Constante temps descente (T2) = ' num2str(T2,'%5.2f \t(')  num2str(T2_ic(1),'%5.2f ; ')  num2str(T2_ic(2),'%5.2f)')], 'FontSize', 14) 

text(0.7, 0.9, ['r² = ' num2str(r2,'%5.4f \t') ], 'FontSize', 16) 

% pause
waitfor(msgbox('NEXT','APPUYER POUR CONTINUER'))
end
% keyboard

FILE = fopen(filename4,'w');

fprintf(FILE,'Ajustement linéaire de la pente\n');
fprintf(FILE,'sur %d points jusqu''au point précédant le maximum\n',DelaiPente);
fprintf(FILE,'pente;r²\n');

h=[];
aaa=[];
H=figure;
for p=1:NP
    indp=B(p):B(p+1);
    tp=t(indp);
    yp=y(indp);
    t1=R(p).t1;
    it1=find(tp>=t1,1,'first');
    it0=max(it1-DelaiPente-10,1);
    range=(it0:it1);
    itm=find(yp(range)==max(yp(range)),1,'last')+it0-1;
    itd=max(1,itm-DelaiPente);
    rangefit=itd:(itm-1);
    X=tp(rangefit);
    Y=yp(rangefit);
    Param=polyfit(X,Y,1);
    Yfit = polyval(Param,X);
    yyy  = polyval(Param,tp);
    yresid = Y - Yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length(Y)-1) * var(Y);
    r2 = 1 - SSresid/SStotal;

    close(h);
    h=figure;plot(tp(range),yp(range),'o')
    hold on;plot(tp(it1),yp(it1),'or','MarkerFaceColor','r')
    plot(tp(itm),yp(itm),'^k','MarkerFaceColor','k')
    plot(X,Yfit,'r')
    fprintf(FILE,'%f;%f\n',Param(1),r2);
    if r2>0.98
        aaa(end+1)=Param(1);
    end
    pause(1)
    

    figure(H);
    hold on
    plot(tp(range)-tp(itm),yp(range)-yp(itm-1),'r')
    plot(tp(range)-tp(itm),yyy(range)-yp(itm-1),'b')
end
aaa=sort(aaa);
L=length(aaa);
B1=max(1,floor(L/4));
B2=min(L,floor(3*L/4));
a_med=median(aaa);
a_mean=mean(aaa(B1:B2));
a_std=std(aaa(B1:B2));
fprintf(FILE,'Après élimination des valeurs pour lesquelles r²<0.98\n');
fprintf(FILE,'pente médiane;%f\n',a_med);
fprintf(FILE,'pente moyenne (sur les 2 quartiles médians);%f\n',a_mean);
fprintf(FILE,'écart-type (sur les 2 quartiles médians);%f\n',a_std);


fprintf(FILE,'\n\n%s\n',Tfonc);
fprintf(FILE,'pente a;');
fprintf(FILE,'%f;;',R.a);
fprintf(FILE,'\n');
fprintf(FILE,'a_intervalle de confiance (inf-sup);');
fprintf(FILE,'%f;',R.a_ic);
fprintf(FILE,'\n');
fprintf(FILE,'offset c;');
fprintf(FILE,'%f;;',R.c);
fprintf(FILE,'\n');
fprintf(FILE,'c_intervalle de confiance (inf-sup);');
fprintf(FILE,'%f;',R.c_ic);
fprintf(FILE,'\n');
fprintf(FILE,'amplitude d;');
fprintf(FILE,'%f;;',R.d);
fprintf(FILE,'\n');
fprintf(FILE,'d_intervalle de confiance (inf-sup);');
fprintf(FILE,'%f;',R.d_ic);
fprintf(FILE,'\n');
fprintf(FILE,'temps de montee t0;');
fprintf(FILE,'%f;;',R.t0);
fprintf(FILE,'\n');
fprintf(FILE,'t0_intervalle de confiance (inf-sup);');
fprintf(FILE,'%f;',R.t0_ic);
fprintf(FILE,'\n');
fprintf(FILE,'temps de descente t1;');
fprintf(FILE,'%f;;',R.t1);
fprintf(FILE,'\n');
fprintf(FILE,'t1_intervalle de confiance (inf-sup);');
fprintf(FILE,'%f;',R.t1_ic);
fprintf(FILE,'\n');
fprintf(FILE,'constante de temps de montee T1;');
fprintf(FILE,'%f;;',R.T1);
fprintf(FILE,'\n');
fprintf(FILE,'T1_intervalle de confiance (inf-sup);');
fprintf(FILE,'%f;',R.T1_ic);
fprintf(FILE,'\n');
fprintf(FILE,'constante de temps de descente T2;');
fprintf(FILE,'%f;;',R.T2);
fprintf(FILE,'\n');
fprintf(FILE,'T2_intervalle de confiance (inf-sup);');
fprintf(FILE,'%f;',R.T2_ic);
fprintf(FILE,'\n');
fprintf(FILE,'r²;');
fprintf(FILE,'%f;;',R.r2);
fprintf(FILE,'\n');
fclose(FILE);

% % subplot(3,1,1)
% % hold on
% % plot(t,y_,'-r', 'LineWidth', 2)
% % plot(t,y,'o','MarkerFaceColor','b')
% % grid
% % xlabel('temps (s)')
% % ylabel('Capno')
% % title(filename,'FontSize',16,'FontWeight','bold','Interpreter','none');
% % 
% % subplot(3,1,2)
% % plot(t,r_,'or', 'MarkerFaceColor','r')
% % grid
% % hold on
% % plot(t,zeros(size(t)),'-','LineWidth',3,'Color',[.5 .5 .5])
% % 
% % xlabel('temps (s)')
% % ylabel('résidus')
% % 
% % 
% % saveas(h,filename3,'pdf');









