%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interface Graphique : "capno_scal"

function capno_scal
    % Create a figure and axes
    f = figure('Visible','off');
    %--------------------------------
    % COURBE
    t0=0; t1=70; t2=100; T1=5; T2=5; c=0; d=5; a = 0.001;
    t = [0:1:t2];
    monte =((d-c)-a*(t1-t)).*(t>t0).*(t<=t1).*(1-exp(-(t-t0)/T1));
    descente =(d-c)*exp(-(t-t1)/T2).*(t>t1);
    courbe = monte + descente;
    %--------------------------------
    % SCALOGRAM
    data=courbe';
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
    %--------------------------------
    % GRAPHIQUE
    subplot(2,1,1);
    imagesc(resultat,[-5 5]); caxis([0,2]); colorbar; colormap('jet');
    xlabel('n (temps)'); 
    ylabel('k (echelle)');
    subplot(2,1,2);
    plot(courbe,'r.-');
    h = max(courbe)+max(courbe)*0.1;
    ylim([0 h]); 
    grid on;
    
    
    % Create pop-up menu
    popup = uicontrol('Style', 'popup',...
           'String', {'jet','parula','hsv','hot','cool','gray'},...
           'Position', [20 340 100 50],...
           'Callback', @setmap);    
          
   % Create slider : "a"
    sld = uicontrol('Style', 'slider',...
        'Min',1,'Max',100,'Value',5,...
        'Position', [400 20 120 20],...
        'Callback', @lim_a); 
					
    % Add a text uicontrol to label the slider.
    txt = uicontrol('Style','text',...
        'Position',[400 40 120 20],...
        'String','variable "a"');
    
       % Create slider : "T1"
    sld = uicontrol('Style', 'slider',...
        'Min',1,'Max',40,'Value',5,...
        'Position', [100 20 120 20],...
        'Callback', @lim_T1); 
					
    % Add a text uicontrol to label the slider.
    txt = uicontrol('Style','text',...
        'Position',[100 40 120 20],...
        'String','variable "T1"');
    
    % Make figure visble after adding all components
    f.Visible = 'on';
    % This code uses dot notation to set properties. 
    % Dot notation runs in R2014b and later.
    % For R2014a and earlier: set(f,'Visible','on');
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function setmap(source,event)
        val = source.Value;
        maps = source.String;
        % For R2014a and earlier: 
        % val = get(source,'Value');
        % maps = get(source,'String'); 

        newmap = maps{val};
        colormap(newmap);
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function lim_a(source,event)
        val_a = source.Value/1000;
        % For R2014a and earlier:
        % val = 51 - get(source,'Value');

        a = val_a;
        monte =((d-c)-a*(t1-t)).*(t>t0).*(t<=t1).*(1-exp(-(t-t0)/T1));
        descente =(d-c)*exp(-(t-t1)/T2).*(t>t1);
        courbe = monte + descente;
        %--------------------------------
        % SCALOGRAM
        data=courbe';
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
        %--------------------------------
        % GRAPHIQUE
        subplot(2,1,1);
        imagesc(resultat,[-5 5]); caxis([0,2]); colorbar; colormap('jet');
        xlabel('n (temps)'); 
        ylabel('k (echelle)');
        subplot(2,1,2);
        plot(courbe,'r.-');
        h = max(courbe)+max(courbe)*0.1;
        ylim([0 h]); 
        grid on;
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function lim_T1(source,event)
        val_T1 = source.Value;
        % For R2014a and earlier:
        % val = 51 - get(source,'Value');

        T1 = val_T1;
        monte =((d-c)-a*(t1-t)).*(t>t0).*(t<=t1).*(1-exp(-(t-t0)/T1));
        descente =(d-c)*exp(-(t-t1)/T2).*(t>t1);
        courbe = monte + descente;
        %--------------------------------
        % SCALOGRAM
        data=courbe';
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
        %--------------------------------
        % GRAPHIQUE
        subplot(2,1,1);
        imagesc(resultat,[-5 5]); caxis([0,2]); colorbar; colormap('jet');
        xlabel('n (temps)'); 
        ylabel('k (echelle)');
        subplot(2,1,2);
        plot(courbe,'r.-');
        h = max(courbe)+max(courbe)*0.1;
        ylim([0 h]); 
        grid on;
    end

end



