%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Travail "Transformée en Ondelette et Capnographie"
% Sébastien Harscoat le 19/06/2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%% type de fonction :

% load_CAPNO();
% scalogram_CAPNO(y);
% imagesc(bpco,[-5 5]); caxis([0,2]); colorbar;   % affichage du scalogram
% figure_SCALOGRAM(courbe,resultat);      % affiche scalogram et courbe
% global_CAPNO(t1,t2, alpha, beta, a);    % produit capnogram et scalogram
% diff_SCALOGRAM(courbe1,courbe2,resultat1,resultat2,fig='on');




%% initialisation
clc
clear
close all



%% production de courbes et leur scalogram

% pour rappel le programme "global_CAPNO"
% donne en sortie une courbe et son scalogram
% et prend en entrée les paramètres suivant : t1, t2, alpha, beta, a
% function[courbe,resultat] = global_CAPNO(t1,t2, alpha, beta, a)


%% Production de capnogram avec différente valeur de t1
figure(1);
[courbe_03,resultat_03] = global_CAPNO(70,100,0.03,0.5,0.8);
figure_SCALOGRAM(courbe_03,resultat_03);

figure(2);
[courbe_08,resultat_08] = global_CAPNO(70,100,0.08,0.5,0.8);
figure_SCALOGRAM(courbe_08,resultat_08);

figure(3);
[courbe_15,resultat_15] = global_CAPNO(70,100,0.15,0.5,0.8);
figure_SCALOGRAM(courbe_15,resultat_15);

figure(4);
[courbe_30,resultat_30] = global_CAPNO(70,100,0.30,0.5,0.8);
figure_SCALOGRAM(courbe_30,resultat_30);


% différence entre scalogram avec t1 de référence à 0.30 et les autres
% (0.15, 0.08, 0.03)

diff_result_30_15 = diff_SCALOGRAM(courbe_30,courbe_15,resultat_30,resultat_15);

diff_result_30_08 = diff_SCALOGRAM(courbe_30,courbe_08,resultat_30,resultat_08);

diff_result_30_03 = diff_SCALOGRAM(courbe_30,courbe_03,resultat_30,resultat_03);


% différence entre scalogram avec t1 de référence à 0.15 et les autres
% (0.30, 0.08, 0.03)

diff_result_15_30 = diff_SCALOGRAM(courbe_15,courbe_30,resultat_15,resultat_30);

diff_result_15_08 = diff_SCALOGRAM(courbe_15,courbe_08,resultat_15,resultat_08);

diff_result_15_03 = diff_SCALOGRAM(courbe_15,courbe_03,resultat_15,resultat_03);


% différence entre scalogram avec t1 de référence à 0.08 et les autres
% (0.03)

diff_result_08_03 = diff_SCALOGRAM(courbe_08,courbe_03,resultat_08,resultat_03);


%% Courbe de capnographie de patient "BPCO"


% chargement courbe de capnographie
yA = load_CAPNO();
yA = yA(12:end);

% prodution du scalogram
bpco = scalogram_CAPNO(yA');

% affichage du scalogram et de la courbe
figure();
figure_SCALOGRAM(yA,bpco);

x = 0:length(yA)-1;
xx = 0:0.25:length(yA)-1;
yAA = spline(x,yA,xx);


% prodution du scalogram
bpco_2 = scalogram_CAPNO(yAA);

% affichage du scalogram et de la courbe
figure();
figure_SCALOGRAM(yAA,bpco_2);


figure();
subplot(1,2,1);
figure_SCALOGRAM(yA,bpco);
subplot(1,2,2);
figure_SCALOGRAM(yAA,bpco_2);



%% Courbe de capnographie de patient "BPCO & Emphysémateux"


% chargement courbe de capnographie
yB = load_CAPNO();
yB = yB(15:end);

% prodution du scalogram
bpco_emph = scalogram_CAPNO(yB');

% affichage du scalogram et de la courbe
figure(7);
figure_SCALOGRAM(yB,bpco_emph);






%% Courbe de capnographie de patient à "Identifier"


% chargement courbe de capnographie
yD = load_CAPNO();
yD = yD(8:end);

% prodution du scalogram
Identifier = scalogram_CAPNO(yD');

% affichage du scalogram et de la courbe
figure(8);
figure_SCALOGRAM(yD,Identifier);








%%  Essaie d'animation de courbe

close all

%%%%%%% Construction de courbe avec modification de paramètres %%%%%%%
subplot(1,2,1);
for T1 = 0.5:3:30
    A = courbe_CAPNO(70,100,T1,5,0.003,5);
    %figure; 
    plot(A,'b.-'); ylim([0 6]); grid on;    
    hold on;
    
end

%hold on;

subplot(1,2,2);
for a = 0:0.01:0.1 
    A = courbe_CAPNO(70,100,5,5,a,5);
    %figure; 
    plot(A,'r.-'); ylim([0 6]); grid on; 
    hold on;
    
end




%%

%%%%%%% Construction de courbe avec modification du paramètre "T1" %%%%%%%
% de 0.5 à 30 avec pat de 3


figure;
for T1 = 0.5:3:30
    A = courbe_CAPNO(70,100,T1,5,0.003,5);
    %figure; 
    plot(A,'b.-'); ylim([0 6]); grid on;    
    hold on;
    
end

%%%%%%% data_base de courbe et de scalogram %%%%%%%

R_T1 = [];
n2_T1 = [];
n3_T1 = [];
n4_T1 = [];
n5_T1 = [];

for T1 = 1:15  %p=1:NP
    y = courbe_CAPNO(70,100,T1,5,0.003,5);
    R_T1(T1).y = y;
    
    % SCALOGRAM
    data=y';
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
    resultat=[fliplr(scalogramme)]';
    R_T1(T1).scal = resultat;
    %R(T1).n_2 = max(resultat(2,:));
    %R(T1).n_3 = max(resultat(3,:));
    %R(T1).n_4 = max(resultat(4,:));
    %R(T1).n_5 = max(resultat(5,:));
    n2_T1(T1) = max(resultat(2,:));
    n3_T1(T1) = max(resultat(3,:));
    n4_T1(T1) = max(resultat(4,:));
    n5_T1(T1) = max(resultat(5,:));
    
end


%%%%%  évolution des valeurs max des niveau en fonction de T1
figure;
hold on;
plot(n2_T1,'m.-'); plot(n3_T1,'g.-'); ;plot(n4_T1,'r.-'); plot(n5_T1,'b.-');



%%%%%%% Représentation de SCALOGRAM d'une série de courbe %%%%%%%



close all;


for T1 = 1:15
    %figure();
    figure_SCALOGRAM(R_T1(T1).y,R_T1(T1).scal);
   
end





%%

%%%%%%% Construction de courbe avec modification du paramètre "a" %%%%%%%
% de 0 à 0.1 avec pat de 0.01


figure;
for a = 0:0.01:0.1 
    A = courbe_CAPNO(70,100,5,5,a,5);
    %figure; 
    plot(A,'r.-'); ylim([0 6]); grid on; 
    hold on;
    
end

%%%%%%% data_base de courbe et de scalogram %%%%%%%

R_a = [];
n2_a = [];
n3_a = [];
n4_a = [];
n5_a = [];
a = 0:0.01:0.1;

for p = 1:length(a)
    y = courbe_CAPNO(70,100,5,5,a(p),5);
    R_a(p).y = y;
    
    % SCALOGRAM
    data=y';
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
    resultat=[fliplr(scalogramme)]';
    R_a(p).scal = resultat;
    %R(T1).n_2 = max(resultat(2,:));
    %R(T1).n_3 = max(resultat(3,:));
    %R(T1).n_4 = max(resultat(4,:));
    %R(T1).n_5 = max(resultat(5,:));
    n2_a(p) = max(resultat(2,:));
    n3_a(p) = max(resultat(3,:));
    n4_a(p) = max(resultat(4,:));
    n5_a(p) = max(resultat(5,:));
    
end


%%%%%  évolution des valeurs max des niveau en fonction de T1
figure;
hold on;
plot(n2_a,'m.-'); plot(n3_a,'g.-'); ;plot(n4_a,'r.-'); plot(n5_a,'b.-');



%%%%%%% Représentation de SCALOGRAM d'une série de courbe %%%%%%%



close all;

for p = 1:length(a)
    %figure();
    figure_SCALOGRAM(R_a(p).y,R_a(p).scal);
    
end





%%

%%%%%%% Construction de courbe avec modification du paramètre "d" %%%%%%%
% de 0 à 10 avec pat de 1


figure;
for d = 0:1:10 
    A = courbe_CAPNO(70,100,5,5,0.001,d);
    %figure; 
    plot(A,'r.-'); ylim([0 11]); grid on; 
    hold on;
    
end

%%%%%%% data_base de courbe et de scalogram %%%%%%%

R_d = [];
n2_d = [];
n3_d = [];
n4_d = [];
n5_d = [];
d = 0:1:10;

for p = 1:length(d)
    y = courbe_CAPNO(70,100,5,5,0.001,d(p));
    R_d(p).y = y;
    
    % SCALOGRAM
    data=y';
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
    resultat=[fliplr(scalogramme)]';
    R_d(p).scal = resultat;
    %R(T1).n_2 = max(resultat(2,:));
    %R(T1).n_3 = max(resultat(3,:));
    %R(T1).n_4 = max(resultat(4,:));
    %R(T1).n_5 = max(resultat(5,:));
    n2_d(p) = max(resultat(2,:));
    n3_d(p) = max(resultat(3,:));
    n4_d(p) = max(resultat(4,:));
    n5_d(p) = max(resultat(5,:));
    
end


%%%%%  évolution des valeurs max des niveau en fonction de T1
figure;
hold on;
plot(n2_d,'m.-'); plot(n3_d,'g.-'); ;plot(n4_d,'r.-'); plot(n5_d,'b.-');



%%%%%%% Représentation de SCALOGRAM d'une série de courbe %%%%%%%



close all;

for p = 1:length(d)
    %figure();
    %figure_SCALOGRAM(R_d(p).y,R_d(p).scal);
    figure();
    %title('fichier :');
    subplot(2,1,1);
    imagesc(R_d(p).scal,[-5 5]); caxis([0,2]); colormap('jet'); colorbar;
    xlabel('n (temps)'); 
    ylabel('k (echelle)');

    subplot(2,1,2);
    plot(R_d(p).y,'r.-');
    ylim([0 11]); 
    grid on;
    
end



%%

%%%%%%% Construction de courbe avec modification du paramètre "T2" %%%%%%%
% de 0.5 à 30 avec pat de 3


figure;
for T2 = 0.5:3:30
    A = courbe_CAPNO(70,100,5,T2,0.003,5);
    %figure; 
    plot(A,'b.-'); ylim([0 6]); grid on;    
    hold on;
    
end

%%%%%%% data_base de courbe et de scalogram %%%%%%%

R_T2 = [];
n2_T2 = [];
n3_T2 = [];
n4_T2 = [];
n5_T2 = [];

for T2 = 1:15  %p=1:NP
    y = courbe_CAPNO(70,100,5,T2,0.003,5);
    R_T2(T2).y = y;
    
    % SCALOGRAM
    data=y';
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
    resultat=[fliplr(scalogramme)]';
    R_T1(T2).scal = resultat;
    %R(T1).n_2 = max(resultat(2,:));
    %R(T1).n_3 = max(resultat(3,:));
    %R(T1).n_4 = max(resultat(4,:));
    %R(T1).n_5 = max(resultat(5,:));
    n2_T2(T2) = max(resultat(2,:));
    n3_T2(T2) = max(resultat(3,:));
    n4_T2(T2) = max(resultat(4,:));
    n5_T2(T2) = max(resultat(5,:));
    
end


%%%%%  évolution des valeurs max des niveau en fonction de T1
figure;
hold on;
plot(n2_T2,'m.-'); plot(n3_T2,'g.-'); ;plot(n4_T2,'r.-'); plot(n5_T2,'b.-');



%%%%%%% Représentation de SCALOGRAM d'une série de courbe %%%%%%%



close all;


for T2 = 1:15
    %figure();
    figure_SCALOGRAM(R_T2(T2).y,R_T1(T2).scal);
   
end







%%

%%%%%%% Construction de courbe avec modification du paramètre "t1" %%%%%%%
% de 10 à 90 avec pat de 10


figure;
for t1 = 10:10:90 
    A = courbe_CAPNO(t1,100,5,5,0.003,5);
    %figure; 
    plot(A,'r.-'); ylim([0 7]); grid on; 
    hold on;
    
end

%%%%%%% data_base de courbe et de scalogram %%%%%%%

R_t1 = [];
n2_t1 = [];
n3_t1 = [];
n4_t1 = [];
n5_t1 = [];
t1 = 10:10:90;

for p = 1:length(t1)
    y = courbe_CAPNO(t1(p),100,5,5,0.003,5);
    R_t1(p).y = y;
    
    % SCALOGRAM
    data=y';
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
    resultat=[fliplr(scalogramme)]';
    R_t1(p).scal = resultat;
    %R(T1).n_2 = max(resultat(2,:));
    %R(T1).n_3 = max(resultat(3,:));
    %R(T1).n_4 = max(resultat(4,:));
    %R(T1).n_5 = max(resultat(5,:));
    n2_t1(p) = max(resultat(2,:));
    n3_t1(p) = max(resultat(3,:));
    n4_t1(p) = max(resultat(4,:));
    n5_t1(p) = max(resultat(5,:));
    
end


%%%%%  évolution des valeurs max des niveau en fonction de T1
figure;
hold on;
plot(n2_t1,'m.-'); plot(n3_t1,'g.-'); ;plot(n4_t1,'r.-'); plot(n5_t1,'b.-');



%%%%%%% Représentation de SCALOGRAM d'une série de courbe %%%%%%%



close all;

for p = 1:length(t1)
    figure();
    %title('fichier :');
    subplot(2,1,1);
    imagesc(R_t1(p).scal,[-5 5]); caxis([0,2]); colormap('jet'); colorbar;
    xlabel('n (temps)'); 
    ylabel('k (echelle)');

    subplot(2,1,2);
    plot(R_t1(p).y,'r.-');
    ylim([0 7]); 
    grid on;
    
end



%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  courbe de variation de "t1" avec différent paramètre de "a"     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%% Construction de courbe avec modification du paramètre "t1" %%%%%%%
% de 10 à 90 avec pat de 10
% pour une valeur de "a" = 0.003


%%%%%%% data_base de courbe et de scalogram %%%%%%%

R_t1 = [];
n2_t1 = [];
n3_t1 = [];
n4_t1 = [];
n5_t1 = [];
t1 = 10:10:90;

for p = 1:length(t1)
    y = courbe_CAPNO(t1(p),100,5,5,0.003,5);
    R_t1(p).y = y;
    
    % SCALOGRAM
    data=y';
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
    resultat=[fliplr(scalogramme)]';
    R_t1(p).scal = resultat;
    %R(T1).n_2 = max(resultat(2,:));
    %R(T1).n_3 = max(resultat(3,:));
    %R(T1).n_4 = max(resultat(4,:));
    %R(T1).n_5 = max(resultat(5,:));
    n2_t1(p) = max(resultat(2,:));
    n3_t1(p) = max(resultat(3,:));
    n4_t1(p) = max(resultat(4,:));
    n5_t1(p) = max(resultat(5,:));
    
end


%%%%%  évolution des valeurs max des niveau en fonction de T1
figure;
hold on;
plot(n2_t1,'m.-'); plot(n3_t1,'g.-'); ;plot(n4_t1,'r.-'); plot(n5_t1,'b.-');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% Construction de courbe avec modification du paramètre "t1" %%%%%%%
% de 10 à 90 avec pat de 10
% pour une valeur de "a" = 0.015


%%%%%%% data_base de courbe et de scalogram %%%%%%%

R_t1 = [];
n2_t1 = [];
n3_t1 = [];
n4_t1 = [];
n5_t1 = [];
t1 = 10:10:90;

for p = 1:length(t1)
    y = courbe_CAPNO(t1(p),100,5,5,0.015,5);
    R_t1(p).y = y;
    
    % SCALOGRAM
    data=y';
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
    resultat=[fliplr(scalogramme)]';
    R_t1(p).scal = resultat;
    %R(T1).n_2 = max(resultat(2,:));
    %R(T1).n_3 = max(resultat(3,:));
    %R(T1).n_4 = max(resultat(4,:));
    %R(T1).n_5 = max(resultat(5,:));
    n2_t1(p) = max(resultat(2,:));
    n3_t1(p) = max(resultat(3,:));
    n4_t1(p) = max(resultat(4,:));
    n5_t1(p) = max(resultat(5,:));
    
end


%%%%%  évolution des valeurs max des niveau en fonction de T1
figure;
hold on;
plot(n2_t1,'m.-'); plot(n3_t1,'g.-'); ;plot(n4_t1,'r.-'); plot(n5_t1,'b.-');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% Construction de courbe avec modification du paramètre "t1" %%%%%%%
% de 10 à 90 avec pat de 10
% pour une valeur de "a" = 0.03


%%%%%%% data_base de courbe et de scalogram %%%%%%%

R_t1 = [];
n2_t1 = [];
n3_t1 = [];
n4_t1 = [];
n5_t1 = [];
t1 = 10:10:90;

for p = 1:length(t1)
    y = courbe_CAPNO(t1(p),100,5,5,0.03,5);
    R_t1(p).y = y;
    
    % SCALOGRAM
    data=y';
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
    resultat=[fliplr(scalogramme)]';
    R_t1(p).scal = resultat;
    %R(T1).n_2 = max(resultat(2,:));
    %R(T1).n_3 = max(resultat(3,:));
    %R(T1).n_4 = max(resultat(4,:));
    %R(T1).n_5 = max(resultat(5,:));
    n2_t1(p) = max(resultat(2,:));
    n3_t1(p) = max(resultat(3,:));
    n4_t1(p) = max(resultat(4,:));
    n5_t1(p) = max(resultat(5,:));
    
end


%%%%%  évolution des valeurs max des niveau en fonction de T1
figure;
hold on;
plot(n2_t1,'m.-'); plot(n3_t1,'g.-'); ;plot(n4_t1,'r.-'); plot(n5_t1,'b.-');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% Construction de courbe avec modification du paramètre "t1" %%%%%%%
% de 10 à 90 avec pat de 10
% pour une valeur de "a" = 0.045


%%%%%%% data_base de courbe et de scalogram %%%%%%%

R_t1 = [];
n2_t1 = [];
n3_t1 = [];
n4_t1 = [];
n5_t1 = [];
t1 = 10:10:90;

for p = 1:length(t1)
    y = courbe_CAPNO(t1(p),100,5,5,0.045,5);
    R_t1(p).y = y;
    
    % SCALOGRAM
    data=y';
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
    resultat=[fliplr(scalogramme)]';
    R_t1(p).scal = resultat;
    %R(T1).n_2 = max(resultat(2,:));
    %R(T1).n_3 = max(resultat(3,:));
    %R(T1).n_4 = max(resultat(4,:));
    %R(T1).n_5 = max(resultat(5,:));
    n2_t1(p) = max(resultat(2,:));
    n3_t1(p) = max(resultat(3,:));
    n4_t1(p) = max(resultat(4,:));
    n5_t1(p) = max(resultat(5,:));
    
end


%%%%%  évolution des valeurs max des niveau en fonction de T1
figure;
hold on;
plot(n2_t1,'m.-'); plot(n3_t1,'g.-'); ;plot(n4_t1,'r.-'); plot(n5_t1,'b.-');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% Construction de courbe avec modification du paramètre "t1" %%%%%%%
% de 10 à 90 avec pat de 10
% pour une valeur de "a" = 0.060


%%%%%%% data_base de courbe et de scalogram %%%%%%%

R_t1 = [];
n2_t1 = [];
n3_t1 = [];
n4_t1 = [];
n5_t1 = [];
t1 = 10:10:90;

for p = 1:length(t1)
    y = courbe_CAPNO(t1(p),100,5,5,0.060,5);
    R_t1(p).y = y;
    
    % SCALOGRAM
    data=y';
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
    resultat=[fliplr(scalogramme)]';
    R_t1(p).scal = resultat;
    %R(T1).n_2 = max(resultat(2,:));
    %R(T1).n_3 = max(resultat(3,:));
    %R(T1).n_4 = max(resultat(4,:));
    %R(T1).n_5 = max(resultat(5,:));
    n2_t1(p) = max(resultat(2,:));
    n3_t1(p) = max(resultat(3,:));
    n4_t1(p) = max(resultat(4,:));
    n5_t1(p) = max(resultat(5,:));
    
end


%%%%%  évolution des valeurs max des niveau en fonction de T1
figure;
hold on;
plot(n2_t1,'m.-'); plot(n3_t1,'g.-'); ;plot(n4_t1,'r.-'); plot(n5_t1,'b.-');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% Construction de courbe avec modification du paramètre "t1" %%%%%%%
% de 10 à 90 avec pat de 10
% pour une valeur de "a" = 0.075


%%%%%%% data_base de courbe et de scalogram %%%%%%%

R_t1 = [];
n2_t1 = [];
n3_t1 = [];
n4_t1 = [];
n5_t1 = [];
t1 = 10:10:90;

for p = 1:length(t1)
    y = courbe_CAPNO(t1(p),100,5,5,0.075,5);
    R_t1(p).y = y;
    
    % SCALOGRAM
    data=y';
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
    resultat=[fliplr(scalogramme)]';
    R_t1(p).scal = resultat;
    %R(T1).n_2 = max(resultat(2,:));
    %R(T1).n_3 = max(resultat(3,:));
    %R(T1).n_4 = max(resultat(4,:));
    %R(T1).n_5 = max(resultat(5,:));
    n2_t1(p) = max(resultat(2,:));
    n3_t1(p) = max(resultat(3,:));
    n4_t1(p) = max(resultat(4,:));
    n5_t1(p) = max(resultat(5,:));
    
end


%%%%%  évolution des valeurs max des niveau en fonction de T1
figure;
hold on;
plot(n2_t1,'m.-'); plot(n3_t1,'g.-'); ;plot(n4_t1,'r.-'); plot(n5_t1,'b.-');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% Construction de courbe avec modification du paramètre "t1" %%%%%%%
% de 10 à 90 avec pat de 10
% pour une valeur de "a" = 0.090


%%%%%%% data_base de courbe et de scalogram %%%%%%%

R_t1 = [];
n2_t1 = [];
n3_t1 = [];
n4_t1 = [];
n5_t1 = [];
t1 = 10:10:90;

for p = 1:length(t1)
    y = courbe_CAPNO(t1(p),100,5,5,0.090,5);
    R_t1(p).y = y;
    
    % SCALOGRAM
    data=y';
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
    resultat=[fliplr(scalogramme)]';
    R_t1(p).scal = resultat;
    %R(T1).n_2 = max(resultat(2,:));
    %R(T1).n_3 = max(resultat(3,:));
    %R(T1).n_4 = max(resultat(4,:));
    %R(T1).n_5 = max(resultat(5,:));
    n2_t1(p) = max(resultat(2,:));
    n3_t1(p) = max(resultat(3,:));
    n4_t1(p) = max(resultat(4,:));
    n5_t1(p) = max(resultat(5,:));
    
end


%%%%%  évolution des valeurs max des niveau en fonction de T1
figure;
hold on;
plot(n2_t1,'m.-'); plot(n3_t1,'g.-'); ;plot(n4_t1,'r.-'); plot(n5_t1,'b.-');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% Construction de courbe avec modification du paramètre "t1" %%%%%%%
% de 10 à 90 avec pat de 10
% pour une valeur de "a" = 0.1


%%%%%%% data_base de courbe et de scalogram %%%%%%%

R_t1 = [];
n2_t1 = [];
n3_t1 = [];
n4_t1 = [];
n5_t1 = [];
t1 = 10:10:90;

for p = 1:length(t1)
    y = courbe_CAPNO(t1(p),100,5,5,0.1,5);
    R_t1(p).y = y;
    
    % SCALOGRAM
    data=y';
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
    resultat=[fliplr(scalogramme)]';
    R_t1(p).scal = resultat;
    %R(T1).n_2 = max(resultat(2,:));
    %R(T1).n_3 = max(resultat(3,:));
    %R(T1).n_4 = max(resultat(4,:));
    %R(T1).n_5 = max(resultat(5,:));
    n2_t1(p) = max(resultat(2,:));
    n3_t1(p) = max(resultat(3,:));
    n4_t1(p) = max(resultat(4,:));
    n5_t1(p) = max(resultat(5,:));
    
end


%%%%%  évolution des valeurs max des niveau en fonction de T1
figure;
hold on;
plot(n2_t1,'m.-'); plot(n3_t1,'g.-'); ;plot(n4_t1,'r.-'); plot(n5_t1,'b.-');






%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  courbe de variation de "d" avec différent paramètre de "a"     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% Construction de courbe avec modification du paramètre "d" %%%%%%%
% de 0 à 10 avec pat de 1
% pour une valeur de "a" = 0.003"

%%%%%%% data_base de courbe et de scalogram %%%%%%%

R_d = [];
n2_d = [];
n3_d = [];
n4_d = [];
n5_d = [];
d = 0:1:10;

for p = 1:length(d)
    y = courbe_CAPNO(70,100,5,5,0.003,d(p));
    R_d(p).y = y;
    
    % SCALOGRAM
    data=y';
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
    resultat=[fliplr(scalogramme)]';
    R_d(p).scal = resultat;
    %R(T1).n_2 = max(resultat(2,:));
    %R(T1).n_3 = max(resultat(3,:));
    %R(T1).n_4 = max(resultat(4,:));
    %R(T1).n_5 = max(resultat(5,:));
    n2_d(p) = max(resultat(2,:));
    n3_d(p) = max(resultat(3,:));
    n4_d(p) = max(resultat(4,:));
    n5_d(p) = max(resultat(5,:));
    
end


%%%%%  évolution des valeurs max des niveau en fonction de T1
figure;
hold on;
plot(n2_d,'m.-'); plot(n3_d,'g.-'); ;plot(n4_d,'r.-'); plot(n5_d,'b.-');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% Construction de courbe avec modification du paramètre "d" %%%%%%%
% de 0 à 10 avec pat de 1
% pour une valeur de "a" = 0.015

%%%%%%% data_base de courbe et de scalogram %%%%%%%

R_d = [];
n2_d = [];
n3_d = [];
n4_d = [];
n5_d = [];
d = 0:1:10;

for p = 1:length(d)
    y = courbe_CAPNO(70,100,5,5,0.015,d(p));
    R_d(p).y = y;
    
    % SCALOGRAM
    data=y';
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
    resultat=[fliplr(scalogramme)]';
    R_d(p).scal = resultat;
    %R(T1).n_2 = max(resultat(2,:));
    %R(T1).n_3 = max(resultat(3,:));
    %R(T1).n_4 = max(resultat(4,:));
    %R(T1).n_5 = max(resultat(5,:));
    n2_d(p) = max(resultat(2,:));
    n3_d(p) = max(resultat(3,:));
    n4_d(p) = max(resultat(4,:));
    n5_d(p) = max(resultat(5,:));
    
end


%%%%%  évolution des valeurs max des niveau en fonction de T1
figure;
hold on;
plot(n2_d,'m.-'); plot(n3_d,'g.-'); ;plot(n4_d,'r.-'); plot(n5_d,'b.-');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% Construction de courbe avec modification du paramètre "d" %%%%%%%
% de 0 à 10 avec pat de 1
% pour une valeur de "a" = 0.03

%%%%%%% data_base de courbe et de scalogram %%%%%%%

R_d = [];
n2_d = [];
n3_d = [];
n4_d = [];
n5_d = [];
d = 0:1:10;

for p = 1:length(d)
    y = courbe_CAPNO(70,100,5,5,0.03,d(p));
    R_d(p).y = y;
    
    % SCALOGRAM
    data=y';
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
    resultat=[fliplr(scalogramme)]';
    R_d(p).scal = resultat;
    %R(T1).n_2 = max(resultat(2,:));
    %R(T1).n_3 = max(resultat(3,:));
    %R(T1).n_4 = max(resultat(4,:));
    %R(T1).n_5 = max(resultat(5,:));
    n2_d(p) = max(resultat(2,:));
    n3_d(p) = max(resultat(3,:));
    n4_d(p) = max(resultat(4,:));
    n5_d(p) = max(resultat(5,:));
    
end


%%%%%  évolution des valeurs max des niveau en fonction de T1
figure;
hold on;
plot(n2_d,'m.-'); plot(n3_d,'g.-'); ;plot(n4_d,'r.-'); plot(n5_d,'b.-');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% Construction de courbe avec modification du paramètre "d" %%%%%%%
% de 0 à 10 avec pat de 1
% pour une valeur de "a" = 0.045

%%%%%%% data_base de courbe et de scalogram %%%%%%%

R_d = [];
n2_d = [];
n3_d = [];
n4_d = [];
n5_d = [];
d = 0:1:10;

for p = 1:length(d)
    y = courbe_CAPNO(70,100,5,5,0.045,d(p));
    R_d(p).y = y;
    
    % SCALOGRAM
    data=y';
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
    resultat=[fliplr(scalogramme)]';
    R_d(p).scal = resultat;
    %R(T1).n_2 = max(resultat(2,:));
    %R(T1).n_3 = max(resultat(3,:));
    %R(T1).n_4 = max(resultat(4,:));
    %R(T1).n_5 = max(resultat(5,:));
    n2_d(p) = max(resultat(2,:));
    n3_d(p) = max(resultat(3,:));
    n4_d(p) = max(resultat(4,:));
    n5_d(p) = max(resultat(5,:));
    
end


%%%%%  évolution des valeurs max des niveau en fonction de T1
figure;
hold on;
plot(n2_d,'m.-'); plot(n3_d,'g.-'); ;plot(n4_d,'r.-'); plot(n5_d,'b.-');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% Construction de courbe avec modification du paramètre "d" %%%%%%%
% de 0 à 10 avec pat de 1
% pour une valeur de "a" = 0.06

%%%%%%% data_base de courbe et de scalogram %%%%%%%

R_d = [];
n2_d = [];
n3_d = [];
n4_d = [];
n5_d = [];
d = 0:1:10;

for p = 1:length(d)
    y = courbe_CAPNO(70,100,5,5,0.06,d(p));
    R_d(p).y = y;
    
    % SCALOGRAM
    data=y';
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
    resultat=[fliplr(scalogramme)]';
    R_d(p).scal = resultat;
    %R(T1).n_2 = max(resultat(2,:));
    %R(T1).n_3 = max(resultat(3,:));
    %R(T1).n_4 = max(resultat(4,:));
    %R(T1).n_5 = max(resultat(5,:));
    n2_d(p) = max(resultat(2,:));
    n3_d(p) = max(resultat(3,:));
    n4_d(p) = max(resultat(4,:));
    n5_d(p) = max(resultat(5,:));
    
end


%%%%%  évolution des valeurs max des niveau en fonction de T1
figure;
hold on;
plot(n2_d,'m.-'); plot(n3_d,'g.-'); ;plot(n4_d,'r.-'); plot(n5_d,'b.-');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% Construction de courbe avec modification du paramètre "d" %%%%%%%
% de 0 à 10 avec pat de 1
% pour une valeur de "a" = 0.075

%%%%%%% data_base de courbe et de scalogram %%%%%%%

R_d = [];
n2_d = [];
n3_d = [];
n4_d = [];
n5_d = [];
d = 0:1:10;

for p = 1:length(d)
    y = courbe_CAPNO(70,100,5,5,0.075,d(p));
    R_d(p).y = y;
    
    % SCALOGRAM
    data=y';
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
    resultat=[fliplr(scalogramme)]';
    R_d(p).scal = resultat;
    %R(T1).n_2 = max(resultat(2,:));
    %R(T1).n_3 = max(resultat(3,:));
    %R(T1).n_4 = max(resultat(4,:));
    %R(T1).n_5 = max(resultat(5,:));
    n2_d(p) = max(resultat(2,:));
    n3_d(p) = max(resultat(3,:));
    n4_d(p) = max(resultat(4,:));
    n5_d(p) = max(resultat(5,:));
    
end


%%%%%  évolution des valeurs max des niveau en fonction de T1
figure;
hold on;
plot(n2_d,'m.-'); plot(n3_d,'g.-'); ;plot(n4_d,'r.-'); plot(n5_d,'b.-');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% Construction de courbe avec modification du paramètre "d" %%%%%%%
% de 0 à 10 avec pat de 1
% pour une valeur de "a" = 0.09

%%%%%%% data_base de courbe et de scalogram %%%%%%%

R_d = [];
n2_d = [];
n3_d = [];
n4_d = [];
n5_d = [];
d = 0:1:10;

for p = 1:length(d)
    y = courbe_CAPNO(70,100,5,5,0.09,d(p));
    R_d(p).y = y;
    
    % SCALOGRAM
    data=y';
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
    resultat=[fliplr(scalogramme)]';
    R_d(p).scal = resultat;
    %R(T1).n_2 = max(resultat(2,:));
    %R(T1).n_3 = max(resultat(3,:));
    %R(T1).n_4 = max(resultat(4,:));
    %R(T1).n_5 = max(resultat(5,:));
    n2_d(p) = max(resultat(2,:));
    n3_d(p) = max(resultat(3,:));
    n4_d(p) = max(resultat(4,:));
    n5_d(p) = max(resultat(5,:));
    
end


%%%%%  évolution des valeurs max des niveau en fonction de T1
figure;
hold on;
plot(n2_d,'m.-'); plot(n3_d,'g.-'); ;plot(n4_d,'r.-'); plot(n5_d,'b.-');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% Construction de courbe avec modification du paramètre "d" %%%%%%%
% de 0 à 10 avec pat de 1
% pour une valeur de "a" = 0.1

%%%%%%% data_base de courbe et de scalogram %%%%%%%

R_d = [];
n2_d = [];
n3_d = [];
n4_d = [];
n5_d = [];
d = 0:1:10;

for p = 1:length(d)
    y = courbe_CAPNO(70,100,5,5,0.1,d(p));
    R_d(p).y = y;
    
    % SCALOGRAM
    data=y';
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
    resultat=[fliplr(scalogramme)]';
    R_d(p).scal = resultat;
    %R(T1).n_2 = max(resultat(2,:));
    %R(T1).n_3 = max(resultat(3,:));
    %R(T1).n_4 = max(resultat(4,:));
    %R(T1).n_5 = max(resultat(5,:));
    n2_d(p) = max(resultat(2,:));
    n3_d(p) = max(resultat(3,:));
    n4_d(p) = max(resultat(4,:));
    n5_d(p) = max(resultat(5,:));
    
end


%%%%%  évolution des valeurs max des niveau en fonction de T1
figure;
hold on;
plot(n2_d,'m.-'); plot(n3_d,'g.-'); ;plot(n4_d,'r.-'); plot(n5_d,'b.-');







%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  courbe de variation de "T2" avec différent paramètre de "a"     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%% Construction de courbe avec modification du paramètre "T2" %%%%%%%
% de 1 à 15 avec pat de 1
% pour une valeur de "a" = 0.003


%%%%%%% data_base de courbe et de scalogram %%%%%%%

R_T2 = [];
n2_T2 = [];
n3_T2 = [];
n4_T2 = [];
n5_T2 = [];

for T2 = 1:15  %p=1:NP
    y = courbe_CAPNO(70,100,5,T2,0.003,5);
    R_T2(T2).y = y;
    
    % SCALOGRAM
    data=y';
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
    resultat=[fliplr(scalogramme)]';
    R_T2(T2).scal = resultat;
    %R(T1).n_2 = max(resultat(2,:));
    %R(T1).n_3 = max(resultat(3,:));
    %R(T1).n_4 = max(resultat(4,:));
    %R(T1).n_5 = max(resultat(5,:));
    n2_T2(T2) = max(resultat(2,:));
    n3_T2(T2) = max(resultat(3,:));
    n4_T2(T2) = max(resultat(4,:));
    n5_T2(T2) = max(resultat(5,:));
    
end


%%%%%  évolution des valeurs max des niveau en fonction de T1
figure;
hold on;
plot(n2_T2,'m.-'); plot(n3_T2,'g.-'); ;plot(n4_T2,'r.-'); plot(n5_T2,'b.-');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% Construction de courbe avec modification du paramètre "T2" %%%%%%%
% de 1 à 15 avec pat de 1
% pour une valeur de "a" = 0.015


%%%%%%% data_base de courbe et de scalogram %%%%%%%

R_T2 = [];
n2_T2 = [];
n3_T2 = [];
n4_T2 = [];
n5_T2 = [];

for T2 = 1:15  %p=1:NP
    y = courbe_CAPNO(70,100,5,T2,0.015,5);
    R_T2(T2).y = y;
    
    % SCALOGRAM
    data=y';
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
    resultat=[fliplr(scalogramme)]';
    R_T2(T2).scal = resultat;
    %R(T1).n_2 = max(resultat(2,:));
    %R(T1).n_3 = max(resultat(3,:));
    %R(T1).n_4 = max(resultat(4,:));
    %R(T1).n_5 = max(resultat(5,:));
    n2_T2(T2) = max(resultat(2,:));
    n3_T2(T2) = max(resultat(3,:));
    n4_T2(T2) = max(resultat(4,:));
    n5_T2(T2) = max(resultat(5,:));
    
end


%%%%%  évolution des valeurs max des niveau en fonction de T1
figure;
hold on;
plot(n2_T2,'m.-'); plot(n3_T2,'g.-'); ;plot(n4_T2,'r.-'); plot(n5_T2,'b.-');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% Construction de courbe avec modification du paramètre "T2" %%%%%%%
% de 1 à 15 avec pat de 1
% pour une valeur de "a" = 0.03


%%%%%%% data_base de courbe et de scalogram %%%%%%%

R_T2 = [];
n2_T2 = [];
n3_T2 = [];
n4_T2 = [];
n5_T2 = [];

for T2 = 1:15  %p=1:NP
    y = courbe_CAPNO(70,100,5,T2,0.03,5);
    R_T2(T2).y = y;
    
    % SCALOGRAM
    data=y';
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
    resultat=[fliplr(scalogramme)]';
    R_T2(T2).scal = resultat;
    %R(T1).n_2 = max(resultat(2,:));
    %R(T1).n_3 = max(resultat(3,:));
    %R(T1).n_4 = max(resultat(4,:));
    %R(T1).n_5 = max(resultat(5,:));
    n2_T2(T2) = max(resultat(2,:));
    n3_T2(T2) = max(resultat(3,:));
    n4_T2(T2) = max(resultat(4,:));
    n5_T2(T2) = max(resultat(5,:));
    
end


%%%%%  évolution des valeurs max des niveau en fonction de T1
figure;
hold on;
plot(n2_T2,'m.-'); plot(n3_T2,'g.-'); ;plot(n4_T2,'r.-'); plot(n5_T2,'b.-');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% Construction de courbe avec modification du paramètre "T2" %%%%%%%
% de 1 à 15 avec pat de 1
% pour une valeur de "a" = 0.045


%%%%%%% data_base de courbe et de scalogram %%%%%%%

R_T2 = [];
n2_T2 = [];
n3_T2 = [];
n4_T2 = [];
n5_T2 = [];

for T2 = 1:15  %p=1:NP
    y = courbe_CAPNO(70,100,5,T2,0.045,5);
    R_T2(T2).y = y;
    
    % SCALOGRAM
    data=y';
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
    resultat=[fliplr(scalogramme)]';
    R_T2(T2).scal = resultat;
    %R(T1).n_2 = max(resultat(2,:));
    %R(T1).n_3 = max(resultat(3,:));
    %R(T1).n_4 = max(resultat(4,:));
    %R(T1).n_5 = max(resultat(5,:));
    n2_T2(T2) = max(resultat(2,:));
    n3_T2(T2) = max(resultat(3,:));
    n4_T2(T2) = max(resultat(4,:));
    n5_T2(T2) = max(resultat(5,:));
    
end


%%%%%  évolution des valeurs max des niveau en fonction de T1
figure;
hold on;
plot(n2_T2,'m.-'); plot(n3_T2,'g.-'); ;plot(n4_T2,'r.-'); plot(n5_T2,'b.-');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% Construction de courbe avec modification du paramètre "T2" %%%%%%%
% de 1 à 15 avec pat de 1
% pour une valeur de "a" = 0.06


%%%%%%% data_base de courbe et de scalogram %%%%%%%

R_T2 = [];
n2_T2 = [];
n3_T2 = [];
n4_T2 = [];
n5_T2 = [];

for T2 = 1:15  %p=1:NP
    y = courbe_CAPNO(70,100,5,T2,0.06,5);
    R_T2(T2).y = y;
    
    % SCALOGRAM
    data=y';
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
    resultat=[fliplr(scalogramme)]';
    R_T2(T2).scal = resultat;
    %R(T1).n_2 = max(resultat(2,:));
    %R(T1).n_3 = max(resultat(3,:));
    %R(T1).n_4 = max(resultat(4,:));
    %R(T1).n_5 = max(resultat(5,:));
    n2_T2(T2) = max(resultat(2,:));
    n3_T2(T2) = max(resultat(3,:));
    n4_T2(T2) = max(resultat(4,:));
    n5_T2(T2) = max(resultat(5,:));
    
end


%%%%%  évolution des valeurs max des niveau en fonction de T1
figure;
hold on;
plot(n2_T2,'m.-'); plot(n3_T2,'g.-'); ;plot(n4_T2,'r.-'); plot(n5_T2,'b.-');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% Construction de courbe avec modification du paramètre "T2" %%%%%%%
% de 1 à 15 avec pat de 1
% pour une valeur de "a" = 0.075


%%%%%%% data_base de courbe et de scalogram %%%%%%%

R_T2 = [];
n2_T2 = [];
n3_T2 = [];
n4_T2 = [];
n5_T2 = [];

for T2 = 1:15  %p=1:NP
    y = courbe_CAPNO(70,100,5,T2,0.075,5);
    R_T2(T2).y = y;
    
    % SCALOGRAM
    data=y';
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
    resultat=[fliplr(scalogramme)]';
    R_T2(T2).scal = resultat;
    %R(T1).n_2 = max(resultat(2,:));
    %R(T1).n_3 = max(resultat(3,:));
    %R(T1).n_4 = max(resultat(4,:));
    %R(T1).n_5 = max(resultat(5,:));
    n2_T2(T2) = max(resultat(2,:));
    n3_T2(T2) = max(resultat(3,:));
    n4_T2(T2) = max(resultat(4,:));
    n5_T2(T2) = max(resultat(5,:));
    
end


%%%%%  évolution des valeurs max des niveau en fonction de T1
figure;
hold on;
plot(n2_T2,'m.-'); plot(n3_T2,'g.-'); ;plot(n4_T2,'r.-'); plot(n5_T2,'b.-');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% Construction de courbe avec modification du paramètre "T2" %%%%%%%
% de 1 à 15 avec pat de 1
% pour une valeur de "a" = 0.09


%%%%%%% data_base de courbe et de scalogram %%%%%%%

R_T2 = [];
n2_T2 = [];
n3_T2 = [];
n4_T2 = [];
n5_T2 = [];

for T2 = 1:15  %p=1:NP
    y = courbe_CAPNO(70,100,5,T2,0.09,5);
    R_T2(T2).y = y;
    
    % SCALOGRAM
    data=y';
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
    resultat=[fliplr(scalogramme)]';
    R_T2(T2).scal = resultat;
    %R(T1).n_2 = max(resultat(2,:));
    %R(T1).n_3 = max(resultat(3,:));
    %R(T1).n_4 = max(resultat(4,:));
    %R(T1).n_5 = max(resultat(5,:));
    n2_T2(T2) = max(resultat(2,:));
    n3_T2(T2) = max(resultat(3,:));
    n4_T2(T2) = max(resultat(4,:));
    n5_T2(T2) = max(resultat(5,:));
    
end


%%%%%  évolution des valeurs max des niveau en fonction de T1
figure;
hold on;
plot(n2_T2,'m.-'); plot(n3_T2,'g.-'); ;plot(n4_T2,'r.-'); plot(n5_T2,'b.-');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% Construction de courbe avec modification du paramètre "T2" %%%%%%%
% de 1 à 15 avec pat de 1
% pour une valeur de "a" = 0.1


%%%%%%% data_base de courbe et de scalogram %%%%%%%

R_T2 = [];
n2_T2 = [];
n3_T2 = [];
n4_T2 = [];
n5_T2 = [];

for T2 = 1:15  %p=1:NP
    y = courbe_CAPNO(70,100,5,T2,0.1,5);
    R_T2(T2).y = y;
    
    % SCALOGRAM
    data=y';
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
    resultat=[fliplr(scalogramme)]';
    R_T2(T2).scal = resultat;
    %R(T1).n_2 = max(resultat(2,:));
    %R(T1).n_3 = max(resultat(3,:));
    %R(T1).n_4 = max(resultat(4,:));
    %R(T1).n_5 = max(resultat(5,:));
    n2_T2(T2) = max(resultat(2,:));
    n3_T2(T2) = max(resultat(3,:));
    n4_T2(T2) = max(resultat(4,:));
    n5_T2(T2) = max(resultat(5,:));
    R_T2(T2).evol = [n2_T2; n3_T2; n4_T2; n5_T2];
    
end


%%%%%  évolution des valeurs max des niveau en fonction de T1
figure;
hold on;
plot(n2_T2,'m.-'); plot(n3_T2,'g.-'); ;plot(n4_T2,'r.-'); plot(n5_T2,'b.-');






















