%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Bouton Analyse Temporelle    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function AnalyseTemporelle_pushButton_Callback(hObject, eventdata, handles,varargin)
% hObject    handle to AnalyseTemporelle_pushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
 
editmat=getappdata(0,'editmat');
 
[xorig,fs]=wavread(editmat);
fs
b=getappdata(0,'b')
%b
x=xorig;
length(x)
 
%-------------------------------------------------           Analyse Temporelle du Signal Altere        -----------------------------------------------%
 
 
bankfilter =[1,100,200,300,400,510,630,770,920,1080,1270,1480,1720,2000,2320,2700,3150,3700,4400,5300,6400,7700,9500,12000,15500];
 
 
ms1=fs/1000;                 % maximum speech Fx at 1000Hz
 
ms20=fs/50;                  % minimum speech Fx at 50Hz
 
frfr=bankfilter(1);
frfr
b=getappdata(0,'b')
b
 
comp=1;
if(b>=1)
    while(b>0)
        comp=comp+1;
        comp
        rest=mod(b,10);
        rest
        if(rest==1)
            comp
            loFreq=bankfilter(comp-1);
            loFreq
            hiFreq=bankfilter(comp);
            hiFreq
 
            fB = [loFreq  hiFreq ]/(fs/2);   
            fB
            [c, a] = butter(1, fB, 'stop');
            x = filtfilt(c, a, x);
%            [c,a,cf,ERB,B]=GammaToneMake(fs,16000,loFreq, hiFreq,'moore');                
%            x=filter(c,a,x);
 
        end
 
        b       
        b=(b-rest)/10;
        b
 
        %b=mod((b-rest),10);
 
    end 
 
            wavwrite(x, fs,'stopBand1.wav');        
            %Y=fft(x,length(x));
            %Y1=ifft(abs(Y).*hamming(length(Y)));            
            %t=(0:length(Y1)-1)/fs; 
            %axes(handles.axes_Endommage); 
            %pot(t,real(Y1));     
           % Y1=ifft(abs(x).*hamming(length(x)));
            t=(0:length(x)-1)/fs;
            axes(handles.axes_Endommage);
            %cla;
            set(handles.axes_Endommage,'HandleVisibility','on');
            plot(t,real(x));                    
            xlabel('Time (s)');
            ylabel('Amplitude'); 
 
 
end
 
%---------------------------------           Fin Analyse Temporelle du Signal Altere      -----------------------------------------------%
 
 
%---------------------------------           Analyse Temporelle du Signal Original        -----------------------------------------------%
 
% plot waveform
 
t=(0:length(xorig)-1)/fs;        % times of sampling instants
 
axes(handles.axes_original);
 
%axes_Endommage.cla;
set(handles.axes_Endommage,'HandleVisibility','off');
plot(t,xorig);
 
legend('Waveform');
 
xlabel('Time (s)');
 
ylabel('Amplitude');
 
%---------------------------------           Fin Analyse Temporelle du Signal Original        --------------------------------------------%