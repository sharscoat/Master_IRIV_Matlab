%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Travail sur l'influence de l'échantillonage





%% Courbe de capnographie de patient "BPCO"


% chargement courbe de capnographie
yA = load_CAPNO();
yA = yA(12:end);

% prodution du scalogram
bpco = scalogram_CAPNO(yA');

% affichage du scalogram et de la courbe
figure();
figure_SCALOGRAM(yA,bpco);

x = 1:length(yA);
xx = 1:0.20:length(yA);
yAA = spline(x,yA,xx);    % utilisation de la fonction "spline" de MATALAB


% prodution du scalogram
bpco_2 = scalogram_CAPNO(yAA);

% affichage du scalogram et de la courbe
figure();
figure_SCALOGRAM(yAA,bpco_2);



figure();
figure_SCALOGRAM(yA,bpco);
figure();
figure_SCALOGRAM(yAA,bpco_2);

%figure();
%subplot(1,2,1);
%figure_SCALOGRAM(yA,bpco);
%subplot(1,2,2);
%figure_SCALOGRAM(yAA,bpco_2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure();

%title('N=77  N0=64  n=6      -        N=305  N0=256  n=8');

subplot(2,2,1);
imagesc(bpco); colormap('jet'); colorbar; caxis([0,2.5]);
xlabel('n (temps)'); 
ylabel('k (echelle)');
title('\fontsize{16} N = 77    N0 = 64    n = 6');
subplot(2,2,3);
plot(yA,'r.-');
h = max(yA)+max(yA)*0.1;
ylim([0 h]); 
grid on;
grid minor;

subplot(2,2,2);
imagesc(bpco_2); colormap('jet'); colorbar; caxis([0,2.5]);
xlabel('n (temps)'); 
ylabel('k (echelle)');
title('\fontsize{16} N = 305    N0 = 256    n = 8');
subplot(2,2,4);
plot(yAA,'r.-');
h = max(yAA)+max(yAA)*0.1;
ylim([0 h]); 
grid on;
grid minor;



%title('N=77  N0=64  n=6      -        N=305  N0=256  n=8');
%title('N=77  N0=64  n=6')
%title('N=305  N0=256  n=8')

%% initialisation
clc
clear all
close all


%% constitution de la base de données : "t2" et "t1"

%%%%%%% Construction d'une base de données en matrice des valeurs de Scalogramme
% en fonction des valeurs de "t2" et de "t1" %%%%%%%

%%%%%%% data_base de courbe et de scalogram %%%%%%%

tic

t2 = 50:5:545;

for j=1:100
tt2(j) = t2(j)*(2/3);
t1 = 0:tt2(j)/100:tt2(j);
    for p = 1:100
        y = courbe_CAPNO(t1(p),t2(j),5,5,0.03,5);

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
        
        An2(j,p) = max(resultat(2,:));
        An3(j,p) = max(resultat(3,:));
        An4(j,p) = max(resultat(4,:));
        An5(j,p) = max(resultat(5,:));
    end

end 

toc


%%%%%  évolution des valeurs max des niveau en fonction de "d" et "a"


figure;
hold on;
colormap('jet')
subplot(2,2,1);
mesh(An2); colorbar;colormap('jet');    %zlim([0 1.2]);caxis([0,1.2]);
xlabel('t1'); ylabel('t2');zlabel('z'); title('niveau 2');
subplot(2,2,2);
mesh(An3); colorbar;colormap('jet');    %zlim([0 1.2]);caxis([0,1.2]);
xlabel('t1'); ylabel('t2');zlabel('z'); title('niveau 3');
subplot(2,2,3);
mesh(An4); colorbar;colormap('jet');    %zlim([0 1.2]);caxis([0,1.2]);
xlabel('t1'); ylabel('t2');zlabel('z'); title('niveau 4');
subplot(2,2,4);
mesh(An5); colorbar;colormap('jet');    %zlim([0 1.2]);caxis([0,1.2]);
xlabel('t1'); ylabel('t2');zlabel('z'); title('niveau 5');







%--------------------------------------------------------------------------
%              %%%%%%%%%%% NORMALISATION par N0 %%%%%%%%%%%


%% constitution de la base de données : "t2" et "t1"

%%%%%%% Construction d'une base de données en matrice des valeurs de Scalogramme
% en fonction des valeurs de "t2" et de "t1" %%%%%%%

%%%%%%% data_base de courbe et de scalogram %%%%%%%

tic

t2 = 50:5:545;

for j=1:100
tt2(j) = t2(j)*(2/3);
t1 = 0:tt2(j)/100:tt2(j);
    for p = 1:100
        y = courbe_CAPNO(t1(p),t2(j),5,5,0.003,5);

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
        %resultat=[fliplr(scalogramme)]';
        resultat_1=[fliplr(scalogramme)]';
        [resultat] = resultat_1/sqrt(N0);
        
        An2(j,p) = max(resultat(2,:));
        An3(j,p) = max(resultat(3,:));
        An4(j,p) = max(resultat(4,:));
        An5(j,p) = max(resultat(5,:));
    end

end 

toc


%%%%%  évolution des valeurs max des niveau en fonction de "d" et "a"


figure;
hold on;
colormap('jet')
subplot(2,2,1);
mesh(An2); colorbar;colormap('jet');    %zlim([0 1.2]);caxis([0,1.2]);
xlabel('t1'); ylabel('t2');zlabel('z'); title('niveau 2');
subplot(2,2,2);
mesh(An3); colorbar;colormap('jet');    %zlim([0 1.2]);caxis([0,1.2]);
xlabel('t1'); ylabel('t2');zlabel('z'); title('niveau 3');
subplot(2,2,3);
mesh(An4); colorbar;colormap('jet');    %zlim([0 1.2]);caxis([0,1.2]);
xlabel('t1'); ylabel('t2');zlabel('z'); title('niveau 4');
subplot(2,2,4);
mesh(An5); colorbar;colormap('jet');    %zlim([0 1.2]);caxis([0,1.2]);
xlabel('t1'); ylabel('t2');zlabel('z'); title('niveau 5');







%--------------------------------------------------------------------------
%%  ANALYSE UNIVARIEE 1



%%%%%%% data_base de courbe et de scalogram %%%%%%%

R_t2 = [];
n2_t2 = [];
n3_t2 = [];
n4_t2 = [];
n5_t2 = [];
t2 = [50 100 150 300 700 1300];  

for p = 1:length(t2)
    y = courbe_CAPNO(t2(p)*(2/3),t2(p),5,5,0.003,5);
    x = 1:length(y);
    xx = 1:0.20:length(y);
    yy = spline(x,y,xx);
    R_t2(p).y = y;
    
    % SCALOGRAM
    data=yy';
    data=data/max(abs(data)); 
    N=length(data); 
    
    % conversion en puissance de 2 
    n=floor(log(N)/log(2)); 
    %n=6;
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
    %resultat_1=[fliplr(scalogramme)]';
    %[resultat] = resultat_1/sqrt(N0);
    R_t2(p).scal = resultat;
    n2_t2(p) = max(resultat(2,:));
    n3_t2(p) = max(resultat(3,:));
    n4_t2(p) = max(resultat(4,:));
    n5_t2(p) = max(resultat(5,:));
    
end


%%%%%  évolution des valeurs max des niveau en fonction de T1
figure;
hold on;
plot(n2_t2,'m.-'); plot(n3_t2,'g.-'); ;plot(n4_t2,'r.-'); plot(n5_t2,'b.-')








%--------------------------------------------------------------------------
%%  ANALYSE UNIVARIEE 2



%%%%%%% data_base de courbe et de scalogram %%%%%%%

R_t2 = [];
n2_t2 = [];
n3_t2 = [];
n4_t2 = [];
n5_t2 = []; 

echant = [1 0.75 0.50 0.40 0.30 0.25 0.20 0.10];  % augmentation croissante

for p = 1:length(echant)
    y = courbe_CAPNO(70,100,1,5,0.003,5);
    x = 1:length(y);
    xx = 1:echant(p):length(y);
    yy = spline(x,y,xx);
    R_t2(p).y = yy;
    
    % SCALOGRAM
    data=yy';
    data=data/max(abs(data)); 
    N=length(data); 
    
    % conversion en puissance de 2 
    n=floor(log(N)/log(2)); 
    %n=6;
    %n=n-1;
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
    resultat_1=[fliplr(scalogramme)]';
    [resultat] = resultat_1/sqrt(N0);
    R_t2(p).scal = resultat;
    n2_t2(p) = max(resultat(2,:));
    n3_t2(p) = max(resultat(3,:));
    n4_t2(p) = max(resultat(4,:));
    n5_t2(p) = max(resultat(5,:));
    
end


%%%%%  évolution des valeurs max des niveau en fonction de T1
figure;
hold on;
plot(n2_t2,'m.-'); plot(n3_t2,'g.-'); ;plot(n4_t2,'r.-'); plot(n5_t2,'b.-')


figure;
hold on;
plot(R_t2(1).y,'m.-'); plot(R_t2(2).y,'g.-'); ;plot(R_t2(3).y,'r.-'); plot(R_t2(4).y,'b.-');
plot(R_t2(5).y,'m.-'); plot(R_t2(6).y,'g.-'); ;plot(R_t2(7).y,'r.-'); plot(R_t2(8).y,'b.-');



%--------------------------------------------------------------------------
%%  BROUILLON

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sébastien HARSCOAT le 20/06/2017
% Programme : "courbe_CAPNO"
% permet de produire une courbes de CAPNOGRAM et de réaliser son SCALOGRAM
% associé.



%% Fonction : "courbe_CAPNO" 


function[courbe,resultat] = courbe_CAPNO(t1,t2, T1, T2, a, d)

    t = [0:1:t2];    % t = [0:1:100];
    t0=0;
    c=0;

    % Paramètres
    % t0=0; t1=70; t2=100;
    % T1=10;   % les valeur de "T1" son comprise entre 2 et 20
    % T2=5;
    % a=0.003;  % les valeur de "a" son comprise entre 0 et 0.1
    % c=0; d=5;

    % Equation de la courbe
    courbe_1 = c + ((d-c)-a*(t1-t)).*(t>t0).*(t<=t1).*(1-exp(-(t-t0)/T1)) + (d-c)*exp(-(t-t1)/T2).*(t>t1);
    monte =((d-c)-a*(t1-t)).*(t>t0).*(t<=t1).*(1-exp(-(t-t0)/T1));
    descente =(d-c)*exp(-(t-t1)/T2).*(t>t1);
    courbe = monte + descente;

    
    % SCALOGRAM
    data=courbe';
    data=data/max(abs(data)); 
    N=length(data); 
        %figure(1); 
        %plot(1:N,data,'r.-'); 
        %ylim([0 1.2]); grid on; 

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

end

