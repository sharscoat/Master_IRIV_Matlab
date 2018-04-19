%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interface Graphique : "capno_scal"

function capno_scal
    % Create a figure and axes
    f = figure('Visible','off');
    %--------------------------------
    % COURBE
    t0=0; t2=100; c=0;
    t1=70; T1=5; T2=5; c=0; d=7; a = 0.001;
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
    n2 = max(resultat(2,:));
    n3 = max(resultat(3,:));
    n4 = max(resultat(4,:));
    n5 = max(resultat(5,:));
    %--------------------------------
    % GRAPHIQUE
    newmap = 'jet';
    subplot(2,3,(1:2));
    imagesc(resultat,[-5 5]); caxis([0,2]); colorbar; colormap(newmap);
    xlabel('n (temps)'); 
    ylabel('k (echelle)');
    subplot(2,3,(4:5));
    plot(courbe,'r.-');
    % h = max(courbe)+max(courbe)*0.1; ylim([0 h]);
    ylim([0 10]); xlim([0 100]);
    grid on;
    subplot(2,3,3);
    x = [2 3 4 5];       % x = categorical({'n2','n3','n4','n5'});
    y = [n2 n3 n4 n5];
    h = bar(x,y); ylim([0 1.2]);
    h.FaceColor = [0 0.7 0.5];
    h.EdgeColor = 'w';
    subplot(2,3,6);
    cla;
    axis off;
    text(0, 1.00, ['a = ' num2str(a,'%5.3f \t')], 'FontSize', 14, 'FontWeight', 'bold'); 
    text(0.3, 1.00, ['c = ' num2str(c,'%5.2f \t')], 'FontSize', 14);
    text(0.6, 1.00, ['d = ' num2str(d,'%5.2f \t')], 'FontSize', 14, 'FontWeight', 'bold'); 
    text(0., 0.85, ['T1 = ' num2str(T1,'%5.1f \t')], 'FontSize', 14, 'FontWeight', 'bold'); 
    text(0.3, 0.85, ['T2 = ' num2str(T2,'%5.1f \t')], 'FontSize', 14);
    text(0.3, 0.70, ['t0 = ' num2str(t0,'%5.0f \t')], 'FontSize', 14);
    text(0, 0.70, ['t1 = ' num2str(t1,'%5.0f \t')], 'FontSize', 14, 'FontWeight', 'bold'); 
    text(0.6, 0.70, ['t1 - t0 = ' num2str(t1-t0,'%5.0f \t') ], 'FontSize', 14);
    text(0., 0.50, ['n2 = ' num2str(n2,'%5.2f \t')], 'FontSize', 14);
    text(0., 0.40, ['n3 = ' num2str(n3,'%5.2f \t')], 'FontSize', 14);
    text(0., 0.30, ['n4 = ' num2str(n4,'%5.2f \t')], 'FontSize', 14);
    text(0., 0.20, ['n5 = ' num2str(n5,'%5.2f \t')], 'FontSize', 14);
    
    % Create pop-up menu
    popup = uicontrol('Style', 'popup',...
           'String', {'jet','parula','hsv','hot','cool','gray'},...
           'Position', [20 340 100 50],...
           'Callback', @setmap);    
          
   % Create slider : "a"
    sld = uicontrol('Style', 'slider',...
        'Min',0.001,'Max',100,'Value',5,...
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
    
    % Create slider : "T2"
    sld = uicontrol('Style', 'slider',...
        'Min',1,'Max',40,'Value',5,...
        'Position', [250 20 120 20],...
        'Callback', @lim_T2); 
					
    % Add a text uicontrol to label the slider.
    txt = uicontrol('Style','text',...
        'Position',[250 40 120 20],...
        'String','variable "T2"');
    
    % Create slider : "d"
    sld = uicontrol('Style', 'slider',...
        'Min',3,'Max',10,'Value',7,...
        'Position', [550 20 120 20],...
        'Callback', @lim_d); 
					
    % Add a text uicontrol to label the slider.
    txt = uicontrol('Style','text',...
        'Position',[550 40 120 20],...
        'String','variable "d"');
    
    % Create slider : "t1"
    sld = uicontrol('Style', 'slider',...
        'Min',10,'Max',100,'Value',70,...
        'Position', [700 20 120 20],...
        'Callback', @lim_t1); 
					
    % Add a text uicontrol to label the slider.
    txt = uicontrol('Style','text',...
        'Position',[700 40 120 20],...
        'String','variable "t1"');
    
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
        @newcourbe;
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
        n2 = max(resultat(2,:));
        n3 = max(resultat(3,:));
        n4 = max(resultat(4,:));
        n5 = max(resultat(5,:));
        %--------------------------------
        % GRAPHIQUE
        subplot(2,3,(1:2));
        imagesc(resultat,[-5 5]); caxis([0,2]); colorbar; colormap(newmap);
        xlabel('n (temps)'); 
        ylabel('k (echelle)');
        subplot(2,3,(4:5));
        plot(courbe,'r.-');
        % h = max(courbe)+max(courbe)*0.1; ylim([0 h]);
        ylim([0 10]); xlim([0 100]);
        grid on;
        subplot(2,3,3);
        x = [2 3 4 5];       % x = categorical({'n2','n3','n4','n5'});
        y = [n2 n3 n4 n5];
        h = bar(x,y); ylim([0 1.2]);
        h.FaceColor = [0 0.7 0.5];
        h.EdgeColor = 'w';
        subplot(2,3,6);
        cla;
        axis off;
        text(0, 1.00, ['a = ' num2str(a,'%5.3f \t')], 'FontSize', 14, 'FontWeight', 'bold'); 
        text(0.3, 1.00, ['c = ' num2str(c,'%5.2f \t')], 'FontSize', 14);
        text(0.6, 1.00, ['d = ' num2str(d,'%5.2f \t')], 'FontSize', 14, 'FontWeight', 'bold'); 
        text(0., 0.85, ['T1 = ' num2str(T1,'%5.1f \t')], 'FontSize', 14, 'FontWeight', 'bold'); 
        text(0.3, 0.85, ['T2 = ' num2str(T2,'%5.1f \t')], 'FontSize', 14);
        text(0.3, 0.70, ['t0 = ' num2str(t0,'%5.0f \t')], 'FontSize', 14);
        text(0, 0.70, ['t1 = ' num2str(t1,'%5.0f \t')], 'FontSize', 14, 'FontWeight', 'bold'); 
        text(0.6, 0.70, ['t1 - t0 = ' num2str(t1-t0,'%5.0f \t') ], 'FontSize', 14);
        text(0., 0.50, ['n2 = ' num2str(n2,'%5.2f \t')], 'FontSize', 14);
        text(0., 0.40, ['n3 = ' num2str(n3,'%5.2f \t')], 'FontSize', 14);
        text(0., 0.30, ['n4 = ' num2str(n4,'%5.2f \t')], 'FontSize', 14);
        text(0., 0.20, ['n5 = ' num2str(n5,'%5.2f \t')], 'FontSize', 14);
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function lim_T2(source,event)
        val_T2 = source.Value;
        % For R2014a and earlier:
        % val = 51 - get(source,'Value');

        T2 = val_T2;
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
        n2 = max(resultat(2,:));
        n3 = max(resultat(3,:));
        n4 = max(resultat(4,:));
        n5 = max(resultat(5,:));
        %--------------------------------
        % GRAPHIQUE
        subplot(2,3,(1:2));
        imagesc(resultat,[-5 5]); caxis([0,2]); colorbar; colormap(newmap);
        xlabel('n (temps)'); 
        ylabel('k (echelle)');
        subplot(2,3,(4:5));
        plot(courbe,'r.-');
        % h = max(courbe)+max(courbe)*0.1; ylim([0 h]);
        ylim([0 10]); xlim([0 100]);
        grid on;
        subplot(2,3,3);
        x = [2 3 4 5];       % x = categorical({'n2','n3','n4','n5'});
        y = [n2 n3 n4 n5];
        h = bar(x,y); ylim([0 1.2]);
        h.FaceColor = [0 0.7 0.5];
        h.EdgeColor = 'w';
        subplot(2,3,6);
        cla;
        axis off;
        text(0, 1.00, ['a = ' num2str(a,'%5.3f \t')], 'FontSize', 14, 'FontWeight', 'bold'); 
        text(0.3, 1.00, ['c = ' num2str(c,'%5.2f \t')], 'FontSize', 14);
        text(0.6, 1.00, ['d = ' num2str(d,'%5.2f \t')], 'FontSize', 14, 'FontWeight', 'bold'); 
        text(0., 0.85, ['T1 = ' num2str(T1,'%5.1f \t')], 'FontSize', 14, 'FontWeight', 'bold'); 
        text(0.3, 0.85, ['T2 = ' num2str(T2,'%5.1f \t')], 'FontSize', 14);
        text(0.3, 0.70, ['t0 = ' num2str(t0,'%5.0f \t')], 'FontSize', 14);
        text(0, 0.70, ['t1 = ' num2str(t1,'%5.0f \t')], 'FontSize', 14, 'FontWeight', 'bold'); 
        text(0.6, 0.70, ['t1 - t0 = ' num2str(t1-t0,'%5.0f \t') ], 'FontSize', 14);
        text(0., 0.50, ['n2 = ' num2str(n2,'%5.2f \t')], 'FontSize', 14);
        text(0., 0.40, ['n3 = ' num2str(n3,'%5.2f \t')], 'FontSize', 14);
        text(0., 0.30, ['n4 = ' num2str(n4,'%5.2f \t')], 'FontSize', 14);
        text(0., 0.20, ['n5 = ' num2str(n5,'%5.2f \t')], 'FontSize', 14);
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function lim_d(source,event)
        val_d = source.Value;
        % For R2014a and earlier:
        % val = 51 - get(source,'Value');

        d = val_d;
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
        n2 = max(resultat(2,:));
        n3 = max(resultat(3,:));
        n4 = max(resultat(4,:));
        n5 = max(resultat(5,:));
        %--------------------------------
        % GRAPHIQUE
        subplot(2,3,(1:2));
        imagesc(resultat,[-5 5]); caxis([0,2]); colorbar; colormap(newmap);
        xlabel('n (temps)'); 
        ylabel('k (echelle)');
        subplot(2,3,(4:5));
        plot(courbe,'r.-');
        % h = max(courbe)+max(courbe)*0.1; ylim([0 h]);
        ylim([0 10]); xlim([0 100]);
        grid on;
        subplot(2,3,3);
        x = [2 3 4 5];       % x = categorical({'n2','n3','n4','n5'});
        y = [n2 n3 n4 n5];
        h = bar(x,y); ylim([0 1.2]);
        h.FaceColor = [0 0.7 0.5];
        h.EdgeColor = 'w';
        subplot(2,3,6);
        cla;
        axis off;
        text(0, 1.00, ['a = ' num2str(a,'%5.3f \t')], 'FontSize', 14, 'FontWeight', 'bold'); 
        text(0.3, 1.00, ['c = ' num2str(c,'%5.2f \t')], 'FontSize', 14);
        text(0.6, 1.00, ['d = ' num2str(d,'%5.2f \t')], 'FontSize', 14, 'FontWeight', 'bold'); 
        text(0., 0.85, ['T1 = ' num2str(T1,'%5.1f \t')], 'FontSize', 14, 'FontWeight', 'bold'); 
        text(0.3, 0.85, ['T2 = ' num2str(T2,'%5.1f \t')], 'FontSize', 14);
        text(0.3, 0.70, ['t0 = ' num2str(t0,'%5.0f \t')], 'FontSize', 14);
        text(0, 0.70, ['t1 = ' num2str(t1,'%5.0f \t')], 'FontSize', 14, 'FontWeight', 'bold'); 
        text(0.6, 0.70, ['t1 - t0 = ' num2str(t1-t0,'%5.0f \t') ], 'FontSize', 14);
        text(0., 0.50, ['n2 = ' num2str(n2,'%5.2f \t')], 'FontSize', 14);
        text(0., 0.40, ['n3 = ' num2str(n3,'%5.2f \t')], 'FontSize', 14);
        text(0., 0.30, ['n4 = ' num2str(n4,'%5.2f \t')], 'FontSize', 14);
        text(0., 0.20, ['n5 = ' num2str(n5,'%5.2f \t')], 'FontSize', 14);
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function lim_t1(source,event)
        val_t1 = source.Value;
        % For R2014a and earlier:
        % val = 51 - get(source,'Value');

        t1 = val_t1;
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
        n2 = max(resultat(2,:));
        n3 = max(resultat(3,:));
        n4 = max(resultat(4,:));
        n5 = max(resultat(5,:));
        %--------------------------------
        % GRAPHIQUE
        subplot(2,3,(1:2));
        imagesc(resultat,[-5 5]); caxis([0,2]); colorbar; colormap(newmap);
        xlabel('n (temps)'); 
        ylabel('k (echelle)');
        subplot(2,3,(4:5));
        plot(courbe,'r.-');
        % h = max(courbe)+max(courbe)*0.1; ylim([0 h]);
        ylim([0 10]); xlim([0 100]);
        grid on;
        subplot(2,3,3);
        x = [2 3 4 5];       % x = categorical({'n2','n3','n4','n5'});
        y = [n2 n3 n4 n5];
        h = bar(x,y); ylim([0 1.2]);
        h.FaceColor = [0 0.7 0.5];
        h.EdgeColor = 'w';
        subplot(2,3,6);
        cla;
        axis off;
        text(0, 1.00, ['a = ' num2str(a,'%5.3f \t')], 'FontSize', 14, 'FontWeight', 'bold');
        text(0.3, 1.00, ['c = ' num2str(c,'%5.2f \t')], 'FontSize', 14);
        text(0.6, 1.00, ['d = ' num2str(d,'%5.2f \t')], 'FontSize', 14, 'FontWeight', 'bold'); 
        text(0., 0.85, ['T1 = ' num2str(T1,'%5.1f \t')], 'FontSize', 14, 'FontWeight', 'bold'); 
        text(0.3, 0.85, ['T2 = ' num2str(T2,'%5.1f \t')], 'FontSize', 14);
        text(0.3, 0.70, ['t0 = ' num2str(t0,'%5.0f \t')], 'FontSize', 14);
        text(0, 0.70, ['t1 = ' num2str(t1,'%5.0f \t')], 'FontSize', 14, 'FontWeight', 'bold'); 
        text(0.6, 0.70, ['t1 - t0 = ' num2str(t1-t0,'%5.0f \t') ], 'FontSize', 14);
        text(0., 0.50, ['n2 = ' num2str(n2,'%5.2f \t')], 'FontSize', 14);
        text(0., 0.40, ['n3 = ' num2str(n3,'%5.2f \t')], 'FontSize', 14);
        text(0., 0.30, ['n4 = ' num2str(n4,'%5.2f \t')], 'FontSize', 14);
        text(0., 0.20, ['n5 = ' num2str(n5,'%5.2f \t')], 'FontSize', 14);
    end

    function newcourbe(source,event)
        a = source.Value;
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
        n2 = max(resultat(2,:));
        n3 = max(resultat(3,:));
        n4 = max(resultat(4,:));
        n5 = max(resultat(5,:));
        %--------------------------------
        subplot(2,3,(1:2));
        imagesc(resultat,[-5 5]); caxis([0,2]); colorbar; colormap(newmap);
        xlabel('n (temps)'); 
        ylabel('k (echelle)');
        subplot(2,3,(4:5));
        plot(courbe,'r.-');
        % h = max(courbe)+max(courbe)*0.1; ylim([0 h]);
        ylim([0 10]); xlim([0 100]);
        grid on;
        subplot(2,3,3);
        x = [2 3 4 5];       % x = categorical({'n2','n3','n4','n5'});
        y = [n2 n3 n4 n5];
        h = bar(x,y); ylim([0 1.2]);
        h.FaceColor = [0 0.7 0.5];
        h.EdgeColor = 'w';
        subplot(2,3,6);
        cla;
        axis off;
        text(0, 1.00, ['a = ' num2str(a,'%5.3f \t')], 'FontSize', 14, 'FontWeight', 'bold'); 
        text(0.3, 1.00, ['c = ' num2str(c,'%5.2f \t')], 'FontSize', 14);
        text(0.6, 1.00, ['d = ' num2str(d,'%5.2f \t')], 'FontSize', 14, 'FontWeight', 'bold'); 
        text(0., 0.85, ['T1 = ' num2str(T1,'%5.1f \t')], 'FontSize', 14, 'FontWeight', 'bold'); 
        text(0.3, 0.85, ['T2 = ' num2str(T2,'%5.1f \t')], 'FontSize', 14);
        text(0.3, 0.70, ['t0 = ' num2str(t0,'%5.0f \t')], 'FontSize', 14);
        text(0, 0.70, ['t1 = ' num2str(t1,'%5.0f \t')], 'FontSize', 14, 'FontWeight', 'bold'); 
        text(0.6, 0.70, ['t1 - t0 = ' num2str(t1-t0,'%5.0f \t') ], 'FontSize', 14);
        text(0., 0.50, ['n2 = ' num2str(n2,'%5.2f \t')], 'FontSize', 14);
        text(0., 0.40, ['n3 = ' num2str(n3,'%5.2f \t')], 'FontSize', 14);
        text(0., 0.30, ['n4 = ' num2str(n4,'%5.2f \t')], 'FontSize', 14);
        text(0., 0.20, ['n5 = ' num2str(n5,'%5.2f \t')], 'FontSize', 14);
    end
end



