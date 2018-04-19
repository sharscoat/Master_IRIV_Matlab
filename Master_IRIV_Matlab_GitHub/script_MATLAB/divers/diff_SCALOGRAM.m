%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sébastien HARSCOAT le 15/06/2017
% Programme : "diff_SCALOGRAM"
% permet de faire la différence entre 2 SCALOGRAM de capnogram et
% d'en afficher les résultats



%% Fonction : "diff_SCALOGRAM"


function[diff_resultat] = diff_SCALOGRAM(courbe1,courbe2,resultat1,resultat2)

    diff_resultat = resultat1 - resultat2;
    
    figure();
    title('courbe 1');
    subplot(2,3,1);
    imagesc(resultat1,[-5 5]); caxis([0,2]); colormap('jet'); colorbar;  
    xlabel('n (temps)'); 
    ylabel('k (echelle)');
    subplot(2,3,4);
    plot(courbe1,'r.-'); 
    ylim([0 1.2]); grid on;
    
    title('courbe 2');
    subplot(2,3,2);
    imagesc(resultat2,[-5 5]); caxis([0,2]); colormap('jet'); colorbar; 
    xlabel('n (temps)'); 
    ylabel('k (echelle)');
    subplot(2,3,5);
    plot(courbe2,'r.-'); 
    ylim([0 1.2]); grid on;

    subplot(2,3,3);
    imagesc(diff_resultat,[-5 5]); caxis([-1,1]); colormap('jet'); colorbar;
    subplot(2,3,6);
    mesh(diff_resultat,[-5 5]); caxis([-1,1]); colormap('jet'); colorbar; 

end 
