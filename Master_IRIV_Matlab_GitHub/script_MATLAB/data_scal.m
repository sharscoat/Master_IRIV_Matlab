%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Simulation de courbe avec différente valeur de "a" et de "T1"
% avec recueil des valeurs du Scalogramme
% après transformée en ondelette
% Sébastien Harscoat le 1/07/2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% initialisation
clc
clear all
close all


%% constitution de la base de données : "a" et "T1"

%%%%%%% Construction d'une base de données en matrice des valeurs de Scalogramme
% en fonction des valeurs de "a" et de "T1" %%%%%%%

%%%%%%% data_base de courbe et de scalogram %%%%%%%

tic

a = 0.001:0.001:0.1;
T1 = 0.15:0.15:15;

for j=1:100

    for p = 1:100
        y = courbe_CAPNO(70,100,T1(p),5,a(j),5);

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


%%%%%  évolution des valeurs max des niveau en fonction de "T1" et "a"

figure;
hold on;
colormap('jet')
subplot(2,2,1);
mesh(An2); colorbar;colormap('jet');	%zlim([0 1.2]);caxis([0,1.2]);
xlabel('T1'); ylabel('a');zlabel('z'); title('n 2');
subplot(2,2,2);
mesh(An3); colorbar;colormap('jet');    %zlim([0 1.2]);caxis([0,1.2]);
xlabel('T1'); ylabel('a');zlabel('z'); title('n 3');
subplot(2,2,3);
mesh(An4); colorbar;colormap('jet');	%zlim([0 1.2]);caxis([0,1.2]);
xlabel('T1'); ylabel('a');zlabel('z'); title('n 4');
subplot(2,2,4);
mesh(An5); colorbar;colormap('jet');    %zlim([0 1.2]);caxis([0,1.2]);
xlabel('T1'); ylabel('a');zlabel('z'); title('n 5');






%% constitution de la base de données : "a" et "d"

%%%%%%% Construction d'une base de données en matrice des valeurs de Scalogramme
% en fonction des valeurs de "a" et de "d" %%%%%%%

%%%%%%% data_base de courbe et de scalogram %%%%%%%

tic

a = 0.001:0.001:0.1;
d = 0.10:0.10:10;

for j=1:100

    for p = 1:100
        y = courbe_CAPNO(70,100,5,5,a(j),d(p));

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



%%%%%%%%%%%%%%%%%%% normalisation de la courbe en fonction de "d"

tic

a = 0.001:0.001:0.1;
d = 0.10:0.10:10;

for j=1:100

    for p = 1:100
        y = courbe_CAPNO(70,100,5,5,a(j),d(p));
        % normalisation
        MM = max(y);
        y = y/MM;

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
xlabel('d'); ylabel('a');zlabel('z'); title('niveau 2');
subplot(2,2,2); 
mesh(An3); colorbar;colormap('jet');    %zlim([0 1.2]);caxis([0,1.2]);
xlabel('d'); ylabel('a');zlabel('z'); title('niveau 3');
subplot(2,2,3);
mesh(An4); colorbar;colormap('jet');   %zlim([0 1.2]);caxis([0,1.2]);
xlabel('d'); ylabel('a');zlabel('z'); title('niveau 4');
subplot(2,2,4);
mesh(An5); colorbar;colormap('jet');   %zlim([0 1.2]);caxis([0,1.2]);
xlabel('d'); ylabel('a');zlabel('z'); title('niveau 5');







%% constitution de la base de données : "a" et "t1"

%%%%%%% Construction d'une base de données en matrice des valeurs de Scalogramme
% en fonction des valeurs de "a" et de "t1" %%%%%%%

%%%%%%% data_base de courbe et de scalogram %%%%%%%

tic

a = 0.001:0.001:0.1;
t1 = 1:1:100;

for j=1:100

    for p = 1:100
        y = courbe_CAPNO(t1(p),100,5,5,a(j),5);

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
xlabel('t1'); ylabel('a');zlabel('z'); title('niveau 2');
subplot(2,2,2);
mesh(An3); colorbar;colormap('jet');    %zlim([0 1.2]);caxis([0,1.2]);
xlabel('t1'); ylabel('a');zlabel('z'); title('niveau 3');
subplot(2,2,3);
mesh(An4); colorbar;colormap('jet');    %zlim([0 1.2]);caxis([0,1.2]);
xlabel('t1'); ylabel('a');zlabel('z'); title('niveau 4');
subplot(2,2,4);
mesh(An5); colorbar;colormap('jet');    %zlim([0 1.2]);caxis([0,1.2]);
xlabel('t1'); ylabel('a');zlabel('z'); title('niveau 5');






%% constitution de la base de données : "a" et "T2"

%%%%%%% Construction d'une base de données en matrice des valeurs de Scalogramme
% en fonction des valeurs de "a" et de "T2" %%%%%%%

%%%%%%% data_base de courbe et de scalogram %%%%%%%

tic

a = 0.001:0.001:0.1;
T2 = 0.15:0.15:15;

for j=1:100

    for p = 1:100
        y = courbe_CAPNO(70,100,5,T2(p),a(j),5);  % T2 soit 100 ou 150

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


%%%%%  évolution des valeurs max des niveau en fonction de "T1" et "a"


figure;
hold on;
colormap('jet')
subplot(2,2,1);
mesh(An2); colorbar;colormap('jet');    %zlim([0 3]);caxis([0,3]);
xlabel('T2'); ylabel('a');zlabel('z'); title('niveau 2');
subplot(2,2,2);
mesh(An3); colorbar;colormap('jet');    %zlim([0 1.2]);caxis([0,1.2]);
xlabel('T2'); ylabel('a');zlabel('z'); title('niveau 3');
subplot(2,2,3);
mesh(An4); colorbar;colormap('jet');    %zlim([0 1.2]);caxis([0,1.2]);
xlabel('T2'); ylabel('a');zlabel('z'); title('niveau 4');
subplot(2,2,4);
mesh(An5); colorbar;colormap('jet');    %zlim([0 1.2]);caxis([0,1.2]);
xlabel('T2'); ylabel('a');zlabel('z'); title('niveau 5');





%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% constitution de la base de données : "T1" et "d"

%%%%%%% Construction d'une base de données en matrice des valeurs de Scalogramme
% en fonction des valeurs de "T1" et de "d" %%%%%%%

%%%%%%% data_base de courbe et de scalogram %%%%%%%

tic

T1 = 0:0.5:50;
d = 0.10:0.10:10;

for j=1:100

    for p = 1:100
        y = courbe_CAPNO(70,100,T1(j),5,0.03,d(p));

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
xlabel('d'); ylabel('T1');zlabel('z'); title('niveau 2');
subplot(2,2,2);
mesh(An3); colorbar;colormap('jet');    %zlim([0 1.2]);caxis([0,1.2]);
xlabel('d'); ylabel('T1');zlabel('z'); title('niveau 3');
subplot(2,2,3);
mesh(An4); colorbar;colormap('jet');   %zlim([0 1.2]);caxis([0,1.2]);
xlabel('d'); ylabel('T1');zlabel('z'); title('niveau 4');
subplot(2,2,4);
mesh(An5); colorbar;colormap('jet');   %zlim([0 1.2]);caxis([0,1.2]);
xlabel('d'); ylabel('T1');zlabel('z'); title('niveau 5');





%%%%%%%%%%%%%%%%%%% normalisation de la courbe en fonction de "d"

tic

T1 = 0.15:0.15:15;
d = 0.10:0.10:10;

for j=1:100

    for p = 1:100
        y = courbe_CAPNO(70,100,T1(j),5,0.03,d(p));
        % normalisation
        MM = max(y);
        y = y/MM;

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
        %resultat_1=[fliplr(scalogramme)]';
        %[resultat] = resultat_1/sqrt(N0);
        
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
mesh(An2); colorbar;colormap('jet');   %zlim([0 1.2]);caxis([0,1.2]);
xlabel('d'); ylabel('T1');zlabel('z'); title('niveau 2');
subplot(2,2,2);
mesh(An3); colorbar;colormap('jet');    %zlim([0 1.2]);caxis([0,1.2]);
xlabel('d'); ylabel('T1');zlabel('z'); title('niveau 3');
subplot(2,2,3);
mesh(An4); colorbar;colormap('jet');   %zlim([0 1.2]);caxis([0,1.2]);
xlabel('d'); ylabel('T1');zlabel('z'); title('niveau 4');
subplot(2,2,4);
mesh(An5); colorbar;colormap('jet');   %zlim([0 1.2]);caxis([0,1.2]);
xlabel('d'); ylabel('T1');zlabel('z'); title('niveau 5');





%% constitution de la base de données : "T1" et "t1"

%%%%%%% Construction d'une base de données en matrice des valeurs de Scalogramme
% en fonction des valeurs de "T1" et de "t1" %%%%%%%

%%%%%%% data_base de courbe et de scalogram %%%%%%%

tic

T1 = 0.15:0.15:15;
t1 = 1:1:100;

for j=1:100

    for p = 1:100
        y = courbe_CAPNO(t1(p),100,T1(j),5,0.03,5);

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
xlabel('t1'); ylabel('T1');zlabel('z'); title('niveau 2');
subplot(2,2,2);
mesh(An3); colorbar;colormap('jet');    %zlim([0 1.2]);caxis([0,1.2]);
xlabel('t1'); ylabel('T1');zlabel('z'); title('niveau 3');
subplot(2,2,3);
mesh(An4); colorbar;colormap('jet');    %zlim([0 1.2]);caxis([0,1.2]);
xlabel('t1'); ylabel('T1');zlabel('z'); title('niveau 4');
subplot(2,2,4);
mesh(An5); colorbar;colormap('jet');    %zlim([0 1.2]);caxis([0,1.2]);
xlabel('t1'); ylabel('T1');zlabel('z'); title('niveau 5');






%% constitution de la base de données : "a" et "T2"

%%%%%%% Construction d'une base de données en matrice des valeurs de Scalogramme
% en fonction des valeurs de "a" et de "T2" %%%%%%%

%%%%%%% data_base de courbe et de scalogram %%%%%%%

tic

T1 = 0.15:0.15:15;     %T1 = 0:0.5:50;
T2 = 0.15:0.15:15;

for j=1:100

    for p = 1:100
        y = courbe_CAPNO(70,100,T1(j),T2(p),0.03,5);  % T2 soit 100 ou 150

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


%%%%%  évolution des valeurs max des niveau en fonction de "T1" et "T2"


figure;
hold on;
colormap('jet')
subplot(2,2,1);
mesh(An2); colorbar;colormap('jet');    %zlim([0 3]);caxis([0,3]);
xlabel('T2'); ylabel('T1');zlabel('z'); title('niveau 2');
subplot(2,2,2);
mesh(An3); colorbar;colormap('jet');    %zlim([0 1.2]);caxis([0,1.2]);
xlabel('T2'); ylabel('T1');zlabel('z'); title('niveau 3');
subplot(2,2,3);
mesh(An4); colorbar;colormap('jet');    %zlim([0 1.2]);caxis([0,1.2]);
xlabel('T2'); ylabel('T1');zlabel('z'); title('niveau 4');
subplot(2,2,4);
mesh(An5); colorbar;colormap('jet');    %zlim([0 1.2]);caxis([0,1.2]);
xlabel('T2'); ylabel('T1');zlabel('z'); title('niveau 5');










%%  valeur de t2

% la valeur de t2 est obligatoirement limitée dans un interval :
% la fréquence respiratoire est comprise entre 10 et 30 cycles par minutes
% en considérant un cycle inspiration et expiration

% pour une fréquence de 10/min le cycle est égale à 6 secondes soit pour un
% pat de temps de 0.02 s, un cycle de 6 secondes correspont à 300 pat

% pour une fréquence de 30/min le cycle est égale à 2 secondes soit pour un
% pat de temps de 0.02 s, un cycle de 2 secondes correspont à 100 pat

% pour une fréquence de 40/min le cycle est égale à 1.5 secondes soit pour un
% pat de temps de 0.02 s, un cycle de 1.5 secondes correspont à 75 pat

% on peut considérer que le cycle respiratoire est compris entre 75 et 300
% pat de temps 




