%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Travail d'analyse des courbes de capnogramme et leur Scalogramme
% Sébastien HARSCOAT  le 15/06/2017




%% initialisation
clc
clear
close all



%% production de courbes et leur scalogram

% pour rappel le programme "global_CAPNO"
% donne en sortie une courbe et son scalogram
% et prend en entrée les paramètres suivant : t1, t2, alpha, beta, a
% function[courbe,resultat] = global_CAPNO(t1,t2, alpha, beta, a)



[courbe_03,resultat_03] = global_CAPNO(70,100,0.03,0.5,0.8);
[courbe_08,resultat_08] = global_CAPNO(70,100,0.08,0.5,0.8);
[courbe_15,resultat_15] = global_CAPNO(70,100,0.15,0.5,0.8);
[courbe_30,resultat_30] = global_CAPNO(70,100,0.30,0.5,0.8);



%% Figure résumé

figure(5);
title('fichier :');
subplot(2,2,1);
imagesc(resultat,[-5 5]); colorbar; 
xlabel('n (temps)'); 
ylabel('k (echelle)');

subplot(2,2,3);
plot(courbe,'r.-'); 
ylim([0 1.2]); grid on;

subplot(2,2,2);
mesh(resultat,[-5 5]);colorbar; 



%% Fonction : "diff_SCALOGRAM"


function[diff_resultat] = diff_SCALOGRAM(courbe1,courbe2,resultat1,resultat2)

    diff_resultat = resultat2 - resultat1;
    
    figure();
    title('courbe 1');
    subplot(2,3,1);
    imagesc(resultat1,[-5 5]); colorbar; 
    xlabel('n (temps)'); 
    ylabel('k (echelle)');
    subplot(2,3,4);
    plot(courbe1,'r.-'); 
    ylim([0 1.2]); grid on;
    
    title('courbe 2');
    subplot(2,3,2);
    imagesc(resultat2,[-5 5]); colorbar; 
    xlabel('n (temps)'); 
    ylabel('k (echelle)');
    subplot(2,3,5);
    plot(courbe2,'r.-'); 
    ylim([0 1.2]); grid on;

    subplot(2,3,3);
    imagesc(diff_resultat,[-5 5]); colorbar;
    subplot(2,3,6);
    mesh(resultat,[-5 5]);colorbar; 

end 



%%

diff_SCALOGRAM(courbe_03,courbe_30,resultat_03,resultat_30);



plot(courbe_08,'b-'); grid on; ylim([0 1.2]);
imagesc(resultat_08,[-5 5]); colorbar; xlabel('n (temps)'); ylabel('k (echelle)');
[courbe_08,resultat_08] = global_CAPNO(70,100,0.08,0.5,0.8);
[courbe_15,resultat_15] = global_CAPNO(70,100,0.15,0.5,0.8);
plot(courbe_15,'b-'); grid on; ylim([0 1.2]);

diff_result_15_08 = resultat_15 - resultat_08



diff_result_15_08 = resultat_15 - resultat_08;
imagesc(diff_result_15_08,[-5 5]); colorbar; xlabel('n (temps)'); ylabel('k (echelle)');
mesh(diff_result_15_08,[-5 5]);colorbar; 
diff_result_08_15 = resultat_08 - resultat_15;
figure(2); mesh(diff_result_08_15,[-5 5]);colorbar;
figure(3);plot(courbe_08,'b-'); figure(4); plot(courbe_15,'r-');


imagesc(diff_result_15_08,[-5 5]); colorbar; xlabel('n (temps)'); ylabel('k (echelle)');
figure(5); imagesc(diff_result_15_08,[-5 5]); colorbar; xlabel('n (temps)'); ylabel('k (echelle)');



















