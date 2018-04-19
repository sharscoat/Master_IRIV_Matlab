%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sébastien HARSCOAT le 16/06/2017
% Programme : "figure_SCALOGRAM"
% permet d'afficher le Scalogram et la courbe de capnographie


%% Fonction : "figure_SCALOGRAM"


function[] = figure_SCALOGRAM(courbe,resultat)

    figure();
    %title('fichier :');
    subplot(2,1,1);
    imagesc(resultat); colormap('jet'); colorbar; caxis([0,0.2]);
    xlabel('n (temps)'); 
    ylabel('k (echelle)');

    subplot(2,1,2);
    plot(courbe,'r.-');
    h = max(courbe)+max(courbe)*0.1;
    ylim([0 h]); 
    grid on;
    
 end