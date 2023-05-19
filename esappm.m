classdef esappm < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                matlab.ui.Figure
        MADEWITHMATLABSTUDENTLICENSEBYERIKASALAZARLabel  matlab.ui.control.Label
        Panel                   matlab.ui.container.Panel
        StackButton             matlab.ui.control.Button
        Step4Label_2            matlab.ui.control.Label
        Step7Label              matlab.ui.control.Label
        Step6Label_2            matlab.ui.control.Label
        Step5Label_2            matlab.ui.control.Label
        Step3Label              matlab.ui.control.Label
        Step2Label              matlab.ui.control.Label
        Step1Label              matlab.ui.control.Label
        DateDatePicker          matlab.ui.control.DatePicker
        DateDatePickerLabel     matlab.ui.control.Label
        PrintReportButton       matlab.ui.control.Button
        UploadFilesButton       matlab.ui.control.Button
        FilterSwitch            matlab.ui.control.RockerSwitch
        FilterSwitchLabel       matlab.ui.control.Label
        HCFHzEditField          matlab.ui.control.NumericEditField
        HCFHzEditFieldLabel     matlab.ui.control.Label
        LCFHzEditField          matlab.ui.control.NumericEditField
        LCFHzEditFieldLabel     matlab.ui.control.Label
        ExpGainEditField        matlab.ui.control.NumericEditField
        ExpGainEditFieldLabel   matlab.ui.control.Label
        DelayEditField          matlab.ui.control.NumericEditField
        DelayEditFieldLabel     matlab.ui.control.Label
        RecordNSpinner          matlab.ui.control.Spinner
        RecordNSpinnerLabel     matlab.ui.control.Label
        LmEditField             matlab.ui.control.NumericEditField
        LmEditFieldLabel        matlab.ui.control.Label
        fsHzEditField           matlab.ui.control.NumericEditField
        fsHzEditFieldLabel      matlab.ui.control.Label
        VcmsEditField           matlab.ui.control.NumericEditField
        VcmsEditFieldLabel      matlab.ui.control.Label
        PileNameEditField       matlab.ui.control.EditField
        PileNameEditFieldLabel  matlab.ui.control.Label
        D3                      matlab.ui.control.UIAxes
        R3                      matlab.ui.control.UIAxes
        R1                      matlab.ui.control.UIAxes
        D1                      matlab.ui.control.UIAxes
        Gain                    matlab.ui.control.UIAxes
    end

    
    properties (Access = private)

        % El presente proyecto desarrollado por Erika Salazar en MATLAB Estudiantil tiene 
        % como objetivo cargar y procesar archivos de texto (.txt). Incluye un proceso de 
        % apilado de los datos cargados, aplicando operaciones de filtrado y normalización. 
        % Se generan visualizaciones gráficas de los datos individuales y apilados, como 
        % reflectogramas y densidad espectral. El proyecto también incorpora elementos de 
        % -interfaz de usuario, como botones y tablas. Se destaca que este desarrollo está 
        % protegido por derechos de autor y queda prohibida su comercialización.

        FileName        % Nombre del archivo que se está cargando
        PathName        % Nombre de la ruta
        filterindex     % Índice del filtro
        data            % Variable que almacena en columnas los datos de los .txt
        stack           % Variable que almacena el apilado
        NumFiles        % Número de archivos cargados
        k               % Variable auxiliar para crear el apilado
        m1              % Ventana con el mensaje 'Cargando datos PIT...'
        hw1             % Desconocido, sin descripción
        samplenum       % Número de muestras
        recordnum       % Número de registros
        stackM          % Variable auxiliar para el apilado
        velocity        % Velocidad en el concreto
        fs              % Frecuencia de muestreo
        PileLength      % Longitud del pilote
        Delay           % Retardo
        ExpGain         % Valor de la ganancia
        Order           % Orden del filtro
        LCF             % Frecuencia inferior del filtro
        HCF             % Frecuencia superior del filtro
        MinPeak         % Variable donde se almacena el pico mínimo
        XMin            % Límite mínimo para el eje X
        XMax            % Límite máximo para el eje X
        YMin            % Límite mínimo para el eje Y
        YMax            % Límite máximo para el eje Y
        dt              % Derivada de tiempo
        dx              % Derivada de posición
        t               % Tiempo
        x               % Posición
        b               % Coeficientes del filtro
        a               % Coeficientes del filtro
        stackMF         % Apilado filtrado
        EndLength       % Longitud final
        EndSample       % Muestra final
        EndTime         % Tiempo final
        TimeDelay       % Retardo de tiempo
        XDelay          % Retardo de posición
        K               % Constante K
        GainFunc        % Función de ganancia
        stackMFG        % Desconocido, sin descripción
        stackMFGN       % Apilado filtrado y normalizado
        Vcms            % Desconocido, sin descripción
        figout          % Figura de salida
        pks             % Picos encontrados
        locs            % Ubicaciones de los picos
        N               % Desconocido, sin descripción
        offsetX         % Desplazamiento en el eje X
        offsetY         % Desplazamiento en el eje Y
        date            % Fecha
        PileName        % Nombre del pilote
        Plot            % Gráfica
        NumFile         % Número de archivo
        DataM           % Datos del archivo seleccionado
        p               % Desconocido, sin descripción
        pxx             % Desconocido, sin descripción
        f               % Frecuencia
        pxxc            % Desconocido, sin descripción
        pB              % Desconocido, sin descripción
        pxxB            % Desconocido, sin descripción
        fB              % Frecuencia B
        pxxcB           % Desconocido, sin descripción
        pC              % Desconocido, sin descripción
        pxxC            % Desconocido, sin descripción
        fC              % Frecuencia C
        pxxcC           % Desconocido, sin descripción
        DateString      % Fecha como cadena de caracteres
        ygain           % Ganancia en el eje Y
        xgain           % Ganancia en el eje X
        i               % Índice
        stackMFN        % Apilado filtrado y normalizado
        tipo            % Desconocido, sin descripción
        s               % Desconocido, sin descripción
        FigoutStack     % Figura para crear la tabla
        UITable         % Tabla con los nombres de los registros cargados por el usuario
        o               % Desconocido, sin descripción
        OKButton        % Elemento de botón
        Step4Label      % Elemento de etiqueta Step4Label
        l               % Desconocido, sin descripción
        UITableData     % Variable donde se almacenan los datos de la tabla
        d               % Desconocido, sin descripción
        FileName2       % Variable donde se guardan los nombres de los archivos seleccionados
        emptyCells      % Variable empleada para eliminar celdas vacías
        inpuFilePath    % Ruta de archivo de entrada
    end

    
    methods (Access = public)
        
        
    end    
        
    
    methods (Access = private)
        
        
    end
      

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: UploadFilesButton
        function UploadFilesButtonPushed(app, event)
        % Solicita al usuario que seleccione archivos TXT y almacena la información en las variables FileName, PathName y filterindex
        [app.FileName, app.PathName, app.filterindex] = uigetfile({'*.*; *.txt', 'Archivos TXT (*.*; *.txt)'}, 'Seleccione archivos txt', 'MultiSelect', 'on');
        % Cambia el directorio de trabajo actual al directorio seleccionado
        cd(app.PathName);
        % Almacena el nombre de archivo seleccionado en la variable FileName
        app.FileName = (app.FileName);        
        % Verifica si el nombre de archivo cumple ciertas condiciones
        if length(app.FileName) < 2 || length(strfind(app.FileName, '.TXT')) < 2 || length(strfind(app.FileName, '.txt')) < 2
            % Carga los datos del archivo en la variable data
            app.data = load(app.FileName);
            app.stack = app.data;            
            % Verifica si es necesario transponer la matriz de datos
            if size(app.data,2) > size(app.data,1)
                app.stack = app.stack';
            end
            % Establece el número de archivos y el contador en 1
            app.NumFiles = 1;
            app.k = 1;
        else
            % Establece el contador en 1 y el número de archivos como la longitud de FileName
            app.k = 1; 
            app.NumFiles = length(app.FileName);            
            % Inicializa la matriz de datos con ceros
            app.data = zeros(length(load(app.FileName{app.k})), length(app.FileName));            
            % Carga los datos del primer archivo en la matriz de datos y en la variable stack
            app.data(:,app.k) = load(app.FileName{app.k});
            app.stack = app.data(:,app.k);            
            % Recorre los archivos restantes y realiza la suma acumulativa en la variable stack
            for k = 2:length(app.FileName)
                app.data(:,k) = load(app.FileName{k});
                app.stack = app.stack + app.data(:,k);
            end
        end
            app.samplenum = length(app.stack); % Número de muestras en el stack
            app.recordnum = size(app.data,2); % Número de registros en los datos
            app.NumFile = app.RecordNSpinner.Value; % Número de archivo seleccionado
            app.DataM = app.data(:,app.NumFile); % Datos del archivo seleccionado
            app.date = app.DateDatePicker.Value; % Fecha seleccionada
            app.DateString = datestr(app.date); % Representación en formato de cadena de la fecha seleccionada
            app.stackM = app.stack/size(app.data,2); % Stack promediado por el número de archivos
            app.stackM = app.stackM/max(abs([min(app.stackM) max(app.stackM)])); % Normalización del stack
            app.PileName = app.PileNameEditField.Value; % Nombre de la pila
            app.velocity = app.VcmsEditField.Value; % Velocidad del concreto en m/s
            app.fs = app.fsHzEditField.Value; % Frecuencia de muestreo en Hz
            app.PileLength = app.LmEditField.Value; % Longitud de la pila en metros
            app.Delay = app.DelayEditField.Value; % Retardo en porcentaje
            app.ExpGain = app.ExpGainEditField.Value; % Ganancia exponencial
            app.Order = 1; % Orden del filtro
            app.LCF = app.LCFHzEditField.Value; % Frecuencia de corte baja en Hz
            app.HCF = app.HCFHzEditField.Value; % Frecuencia de corte alta en Hz
            app.date = app.DateDatePicker.Value; % Fecha seleccionada
            app.MinPeak = 0.75; % Altura mínima de los picos
            app.offsetX = 0.1; % Desplazamiento en el eje X
            app.offsetY = 0.1; % Desplazamiento en el eje Y
            app.XMin = 0.0; % Valor mínimo en el eje X
            app.XMax = 1.2*app.PileLength; % Valor máximo en el eje X
            app.YMin = -1.0; % Valor mínimo en el eje Y
            app.YMax = 1.0; % Valor máximo en el eje Y
            app.dt = 1/app.fs; % Intervalo de tiempo entre muestras
            app.dx = app.dt*app.velocity/2; % Intervalo espacial entre muestras
            app.t = 0:app.dt:(length(app.stack)-1)*app.dt; % Vector de tiempo
            app.x = 0:app.dx:(length(app.stack)-1)*app.dx; % Vector espacial
            [app.b,app.a] = butter(app.Order,[app.LCF app.HCF]/(app.fs/2),'bandpass'); % Coeficientes del filtro Butterworth
            app.stackMF = filter(app.b,app.a,app.stackM); % Aplicación del filtro al stack promediado
            app.EndLength = 1.10*app.PileLength; % Longitud del extremo en metros
            app.EndSample = floor((app.EndLength)/app.dx); % Muestra correspondiente al extremo
            app.EndTime = app.EndSample*app.dt; % Tiempo correspondiente al extremo
            app.TimeDelay = (app.Delay/100)*app.dt*app.EndSample; % Retardo de tiempo
            app.XDelay = (app.Delay/100)*app.dx*app.EndSample; % Retardo espacial
            app.K = log(app.ExpGain)/((app.EndSample*app.dt)-app.TimeDelay); % Cálculo de K
            app.GainFunc = ones(length(app.stackMF),1); % Inicialización de GainFunc como un vector de unos
            for i = 1:length(app.stackMF)
                if app.t(i) < app.TimeDelay
                    app.GainFunc(i) = 1; % Ganancia igual a 1 para tiempos menores que TimeDelay
                elseif app.t(i) >= app.TimeDelay && app.t(app.EndSample) >= app.t(i)
                    app.GainFunc(i) = exp(app.K*(app.t(i) - app.TimeDelay)); % Ganancia exponencial creciente
                else
                    app.GainFunc(i) = app.ExpGain*exp(-app.K*(app.t(i) - app.EndTime)); % Ganancia exponencial decreciente
                end
            end
            app.stackMFN = app.stackMF/max(abs(app.stackMF)); % Normalización del stackMF
            [app.pks,app.locs] = findpeaks(abs(app.stackMFN),'MINPEAKHEIGHT',app.MinPeak); % Encontrar picos en el stackMFN
            
            % Gráficas individuales
            plot(app.R1, app.x,app.DataM, '-b');
            title(app.R1, sprintf('Individual Record Reflectogram. N°= %d', app.NumFile));
            xlim(app.R1,[app.XMin, app.XMax]);
            
            [app.pxx,app.f,app.pxxc] = periodogram(app.DataM,rectwin(length(app.DataM)),length(app.DataM),app.fs,...
            'ConfidenceLevel',0.95);
            plot(app.D1, app.f,10*log10(app.pxx),'-', 'Color','b');
            title(app.D1, sprintf('Individual Record Spectral Density. N°= %d', app.NumFile));
            
            % Gráfica de ganancias
            plot(app.Gain, app.x, app.GainFunc, '-b');
            
            % Gráfica del stack
            plot(app.R3, app.x, app.stack, '-b');
            title (app.R3, 'Stacked Record Reflectogram');
            xlim(app.R3,[app.XMin, app.XMax]);
            
            [app.pxxB,app.fB,app.pxxcB] = periodogram(app.stack,rectwin(length(app.stack)),length(app.stack),app.fs);
            plot(app.D3, app.fB,10*log10(app.pxxB),'-', 'Color','b');
            title (app.D3, 'Spectral Density of the Stacked Record');
        end

        % Value changed function: VcmsEditField
        function VcmsEditFieldValueChanged(app, event)
            % Obtiene el número de archivo seleccionado desde el spinner y lo almacena en la variable NumFile
            app.NumFile = app.RecordNSpinner.Value;      
            % Extrae los datos correspondientes al archivo seleccionado y los almacena en la variable DataM
            app.DataM = app.data(:,app.NumFile);            
            % Obtiene la fecha seleccionada y la convierte en una cadena
            app.date = app.DateDatePicker.Value;
            app.DateString = datestr(app.date );            
            % Calcula la señal apilada normalizada y la almacena en la variable stackM
            app.stackM = app.stack/size(app.data,2);
            app.stackM = app.stackM/max(abs([min(app.stackM) max(app.stackM)]));            
            % Obtiene el nombre de la pila desde el campo de edición y lo almacena en la variable PileName
            app.PileName = app.PileNameEditField.Value;            
            % Obtiene la velocidad del hormigón desde el campo de edición y lo almacena en la variable velocity
            app.velocity = app.VcmsEditField.Value;           
            % Obtiene la frecuencia de muestreo desde el campo de edición y lo almacena en la variable fs
            app.fs = app.fsHzEditField.Value;             
            % Obtiene la longitud de la pila desde el campo de edición y lo almacena en la variable PileLength
            app.PileLength = app.LmEditField.Value;             
            % Obtiene el retardo desde el campo de edición y lo almacena en la variable Delay
            app.Delay = app.DelayEditField.Value;            
            % Obtiene la ganancia exponencial desde el campo de edición y la almacena en la variable ExpGain
            app.ExpGain = app.ExpGainEditField.Value;       
            % Establece el orden del filtro Butterworth en 1
            app.Order = 1;            
            % Obtiene la frecuencia de corte inferior del filtro desde el campo de edición y la almacena en la variable LCF
            app.LCF = app.LCFHzEditField.Value;            
            % Obtiene la frecuencia de corte superior del filtro desde el campo de edición y la almacena en la variable HCF
            app.HCF = app.HCFHzEditField.Value;            
            % Obtiene la fecha seleccionada y la almacena en la variable date
            app.date = app.DateDatePicker.Value;            
            % Establece el valor mínimo para los picos en 0.75
            app.MinPeak = 0.75;            
            % Establece los valores de desplazamiento horizontal y vertical en 0.1
            app.offsetX = 0.1;
            app.offsetY = 0.1;
            % Establece los límites del eje x y eje y en las subfiguras
            app.XMin = 0.0;
            app.XMax = 1.2*app.PileLength;
            app.YMin = -1.0;
            app.YMax = 1.0;            
            % Calcula el paso de tiempo y el paso espacial
            app.dt = 1/app.fs;
            app.dx = app.dt * app.velocity / 2;  % Calcula el espaciado entre posiciones utilizando la velocidad y el intervalo de tiempo
            app.t = 0:app.dt:(length(app.stack)-1)*app.dt;  % Crea un vector de tiempo desde 0 hasta el tiempo final del stack en incrementos de dt
            app.x = 0:app.dx:(length(app.stack)-1)*app.dx;  % Crea un vector de posiciones desde 0 hasta la posición final del stack en incrementos de dx
            [app.b, app.a] = butter(app.Order, [app.LCF app.HCF]/(app.fs/2), 'bandpass');  % Diseña un filtro Butterworth de paso de banda
            app.stackMF = filter(app.b, app.a, app.stackM);  % Aplica el filtro al stack de datos
            app.EndLength = 1.10 * app.PileLength;  % Calcula la longitud final del pilote como un 10% más de la longitud original
            app.EndSample = floor(app.EndLength / app.dx);  % Calcula el número de muestras correspondientes a la longitud final del pilote
            app.EndTime = app.EndSample * app.dt;  % Calcula el tiempo correspondiente a la longitud final del pilote
            app.TimeDelay = (app.Delay/100) * app.dt * app.EndSample;  % Calcula el retardo de tiempo en función del porcentaje de retardo, el intervalo de tiempo y el número de muestras de la longitud final del pilote
            app.XDelay = (app.Delay/100) * app.dx * app.EndSample;  % Calcula el retardo de posición en función del porcentaje de retardo, el espaciado entre posiciones y el número de muestras de la longitud final del pilote
            app.K = log(app.ExpGain) / ((app.EndSample * app.dt) - app.TimeDelay);  % Calcula el valor de K para la función de ganancia, utilizando el logaritmo del valor de ganancia y considerando el tiempo total menos el retardo de tiempo
            app.GainFunc = ones(length(app.stackMF), 1);  % Crea un vector de ganancia inicializado con unos, con una longitud igual al stack filtrado
            for i = 1:length(app.stackMF)
                if app.t(i) < app.TimeDelay
                    app.GainFunc(i) = 1;
                elseif app.t(i) >= app.TimeDelay  && app.t(app.EndSample) >= app.t(i)
                    app.GainFunc(i) = exp(app.K*(app.t(i) - app.TimeDelay));
                else app.t(i) >= app.t(app.EndSample)
                    app.GainFunc(i) = app.ExpGain*exp(-app.K*(app.t(i) - app.EndTime));
                end
            end
            app.stackMFN = app.stackMF/max(abs(app.stackMF));
            [app.pks,app.locs]=findpeaks(abs(app.stackMFN),'MINPEAKHEIGHT',app.MinPeak); 
            %Graficas
            %----Individuales----
            plot(app.R1, app.x,app.DataM, '-b');
            title(app.R1, sprintf('Individual Record Reflectogram. N°= %d', app.NumFile));
            xlim(app.R1,[app.XMin, app.XMax]);
            [app.pxx,app.f,app.pxxc] = periodogram(app.DataM,rectwin(length(app.DataM)),length(app.DataM),app.fs,...
            'ConfidenceLevel',0.95);
            plot(app.D1, app.f,10*log10(app.pxx),'-', Color='b');
             title(app.D1, sprintf('Individual Record Spectral Density. N°= %d', app.NumFile));
            %----Ganancias----
            plot(app.Gain, app.x, app.GainFunc, '-b');
            %----Apilado---
            plot(app.R3, app.x, app.stack, '-b');   %Grafica individual       
            [app.pxxB,app.fB,app.pxxcB] = periodogram(app.stack,rectwin(length(app.stack)),length(app.stack),app.fs);
            plot(app.D3, app.fB,10*log10(app.pxxB),'-', Color='b');
                       
        end

        % Value changed function: fsHzEditField
        function fsHzEditFieldValueChanged(app, event)
        %La función fsHzEditFieldValueChanged es un callback que se activa cuando se cambia el valor
        % del campo fsHzEditField. Esta función realiza diversas operaciones y cálculos relacionados
        % con el procesamiento de datos sísmicos de pilotes. Los pasos principales incluyen la selección 
        % de archivos y datos, normalización y escalado del stack de datos, configuración de parámetros 
        % de filtrado y ganancia, cálculo de intervalos de tiempo y posición, aplicación de un filtro 
        % Butterworth al stack de datos, cálculo de parámetros relacionados con la longitud final del 
        % pilote y el retardo de tiempo, y generación de una función de ganancia. Además, se realizan 
        % gráficas de los datos individuales, ganancia y stack, y se calcula el espectro de frecuencia.
            app.NumFile = app.RecordNSpinner.Value;
            app.DataM = app.data(:,app.NumFile);
            app.date = app.DateDatePicker.Value;
            app.DateString = datestr(app.date );
            app.stackM = app.stack/size(app.data,2);
            app.stackM = app.stackM/max(abs([min(app.stackM) max(app.stackM)]));
            app.PileName = app.PileNameEditField.Value;
            app.velocity = app.VcmsEditField.Value;                 %Concrete velocity m/s
            app.fs = app.fsHzEditField.Value;                       %frequency sampling in Hz
            app.PileLength = app.LmEditField.Value;                 %in meters
            app.Delay = app.DelayEditField.Value;
            app.ExpGain = app.ExpGainEditField.Value;
            app.Order = 1;
            app.LCF = app.LCFHzEditField.Value;
            app.HCF = app.HCFHzEditField.Value;
            app.date = app.DateDatePicker.Value;
            app.MinPeak = 0.75;
            app.offsetX = 0.1;
            app.offsetY = 0.1;
            app.XMin = 0.0;
            app.XMax = 1.2*app.PileLength;
            app.YMin = -1.0;
            app.YMax = 1.0;
            app.dt = 1/app.fs;
            app.dx = app.dt*app.velocity/2;
            app.t = 0:app.dt:(length(app.stack)-1)*app.dt;
            app.x = 0:app.dx:(length(app.stack)-1)*app.dx;
            [app.b,app.a] = butter(app.Order,[app.LCF app.HCF]/(app.fs/2),'bandpass');
            app.stackMF = filter(app.b,app.a,app.stackM);
            app.EndLength = 1.10*app.PileLength; %in meters
            app.EndSample = floor((app.EndLength)/app.dx);
            app.EndTime = app.EndSample*app.dt;
            app.TimeDelay = (app.Delay/100)*app.dt*app.EndSample;
            app.XDelay = (app.Delay/100)*app.dx*app.EndSample;
            app.K = log(app.ExpGain)/((app.EndSample*app.dt)-app.TimeDelay);
            app.GainFunc = ones(length(app.stackMF),1);
            for i = 1:length(app.stackMF)
                if app.t(i) < app.TimeDelay
                    app.GainFunc(i) = 1;
                elseif app.t(i) >= app.TimeDelay  && app.t(app.EndSample) >= app.t(i)
                    app.GainFunc(i) = exp(app.K*(app.t(i) - app.TimeDelay));
                else app.t(i) >= app.t(app.EndSample)
                    app.GainFunc(i) = app.ExpGain*exp(-app.K*(app.t(i) - app.EndTime));
                end
            end
            app.stackMFN = app.stackMF/max(abs(app.stackMF));
            [app.pks,app.locs]=findpeaks(abs(app.stackMFN),'MINPEAKHEIGHT',app.MinPeak); 
            %Graficas
            %----Individuales----
            plot(app.R1, app.x,app.DataM, '-b');
            title(app.R1, sprintf('Individual Record Reflectogram. N°= %d', app.NumFile));
            xlim(app.R1,[app.XMin, app.XMax]);
            [app.pxx,app.f,app.pxxc] = periodogram(app.DataM,rectwin(length(app.DataM)),length(app.DataM),app.fs,...
            'ConfidenceLevel',0.95);
            plot(app.D1, app.f,10*log10(app.pxx),'-', Color='b');
             title(app.D1, sprintf('Individual Record Spectral Density. N°= %d', app.NumFile));
            %----Ganancias----
            yline(app.Gain, app.ExpGain, '-', Color='b');   % Ganancia 
            xlim(app.Gain,[app.XMin, app.PileLength]);
            %----Apilado---
            plot(app.R3, app.x, app.stack, '-b');   %Grafica individual       
            [app.pxxB,app.fB,app.pxxcB] = periodogram(app.stack,rectwin(length(app.stack)),length(app.stack),app.fs);
            plot(app.D3, app.fB,10*log10(app.pxxB),'-', Color='b');
        end

        % Value changed function: RecordNSpinner
        function RecordNSpinnerValueChanged(app, event)
            app.NumFile = app.RecordNSpinner.Value;
            app.DataM = app.data(:,app.NumFile);
            app.date = app.DateDatePicker.Value;
            app.DateString = datestr(app.date );
            app.stackM = app.stack/size(app.data,2);
            app.stackM = app.stackM/max(abs([min(app.stackM) max(app.stackM)]));
            app.PileName = app.PileNameEditField.Value;
            app.velocity = app.VcmsEditField.Value;                 %Concrete velocity m/s
            app.fs = app.fsHzEditField.Value;                       %frequency sampling in Hz
            app.PileLength = app.LmEditField.Value;                 %in meters
            app.Delay = app.DelayEditField.Value;
            app.ExpGain = app.ExpGainEditField.Value;
            app.Order = 1;
            app.LCF = app.LCFHzEditField.Value;
            app.HCF = app.HCFHzEditField.Value;
            app.date = app.DateDatePicker.Value;
            app.MinPeak = 0.75;
            app.offsetX = 0.1;
            app.offsetY = 0.1;
            app.XMin = 0.0;
            app.XMax = 1.2*app.PileLength;
            app.YMin = -1.0;
            app.YMax = 1.0;
            app.dt = 1/app.fs;
            app.dx = app.dt*app.velocity/2;
            app.t = 0:app.dt:(length(app.stack)-1)*app.dt;
            app.x = 0:app.dx:(length(app.stack)-1)*app.dx;
            [app.b,app.a] = butter(app.Order,[app.LCF app.HCF]/(app.fs/2),'bandpass');
            app.stackMF = filter(app.b,app.a,app.stackM);
            app.EndLength = 1.10*app.PileLength; %in meters
            app.EndSample = floor((app.EndLength)/app.dx);
            app.EndTime = app.EndSample*app.dt;
            app.TimeDelay = (app.Delay/100)*app.dt*app.EndSample;
            app.XDelay = (app.Delay/100)*app.dx*app.EndSample;
            app.K = log(app.ExpGain)/((app.EndSample*app.dt)-app.TimeDelay);
            app.GainFunc = ones(length(app.stackMF),1);
            for i = 1:length(app.stackMF)
                if app.t(i) < app.TimeDelay
                    app.GainFunc(i) = 1;
                elseif app.t(i) >= app.TimeDelay  && app.t(app.EndSample) >= app.t(i)
                    app.GainFunc(i) = exp(app.K*(app.t(i) - app.TimeDelay));
                else app.t(i) >= app.t(app.EndSample)
                    app.GainFunc(i) = app.ExpGain*exp(-app.K*(app.t(i) - app.EndTime));
                end
            end
            app.stackMFN = app.stackMF/max(abs(app.stackMF));
            [app.pks,app.locs]=findpeaks(abs(app.stackMFN),'MINPEAKHEIGHT',app.MinPeak); 
            %Graficas
            %----Individuales----
            plot(app.R1, app.x,app.DataM, '-b');
            title(app.R1, sprintf('Individual Record Reflectogram. N°= %d', app.NumFile));
            xlim(app.R1,[app.XMin, app.XMax]);
            [app.pxx,app.f,app.pxxc] = periodogram(app.DataM,rectwin(length(app.DataM)),length(app.DataM),app.fs,...
            'ConfidenceLevel',0.95);
            plot(app.D1, app.f,10*log10(app.pxx),'-', Color='b');
             title(app.D1, sprintf('Individual Record Spectral Density. N°= %d', app.NumFile));
            %----Ganancias----
            plot(app.Gain, app.x, app.GainFunc, '-b');
            %----Apilado---
            plot(app.R3, app.x, app.stack, '-b');   %Grafica individual       
            [app.pxxB,app.fB,app.pxxcB] = periodogram(app.stack,rectwin(length(app.stack)),length(app.stack),app.fs);
            plot(app.D3, app.fB,10*log10(app.pxxB),'-', Color='b');
            
        end

        % Value changed function: FilterSwitch
        function FilterSwitchValueChanged(app, event)
            
            app.t = clock;
            app.samplenum = length(app.stack);
            app.recordnum = size(app.data,2);
            app.stackM = app.stack/size(app.data,2);
            app.stackM = app.stackM/max(abs([min(app.stackM) max(app.stackM)]));
            app.PileName = app.PileNameEditField.Value;
            app.velocity = app.VcmsEditField.Value;                 %Concrete velocity m/s
            app.fs = app.fsHzEditField.Value;                       %frequency sampling in Hz
            app.PileLength = app.LmEditField.Value;                 %in meters
            app.Delay = app.DelayEditField.Value;
            app.ExpGain = app.ExpGainEditField.Value;
            app.Order = 1;
            app.LCF = app.LCFHzEditField.Value;
            app.HCF = app.HCFHzEditField.Value;
            app.date = app.DateDatePicker.Value;
            app.MinPeak = 0.75;
            app.offsetX = 0.1;
            app.offsetY = 0.1;
            app.XMin = 0.0;
            app.XMax = 1.2*app.PileLength;
            app.YMin = -1.0;
            app.YMax = 1.0;
            app.dt = 1/app.fs;
            app.dx = app.dt*app.velocity/2;
            app.t = 0:app.dt:(length(app.stack)-1)*app.dt;
            app.x = 0:app.dx:(length(app.stack)-1)*app.dx;
            [app.b,app.a] = butter(app.Order,[app.LCF app.HCF]/(app.fs/2),'bandpass');
            app.stackMF = filter(app.b,app.a,app.stackM);
            app.EndLength = 1.10*app.PileLength; %in meters
            app.EndSample = floor(app.EndLength/app.dx);
            app.EndTime = app.EndSample*app.dt;
            app.TimeDelay = (app.Delay/100)*app.dt*app.EndSample;
            app.XDelay = (app.Delay/100)*app.dx*app.EndSample;
            app.K = log(app.ExpGain)/((app.EndSample*app.dt)-app.TimeDelay);
            app.stackMFN = app.stackMF/max(abs(app.stackMF));
            [app.pks,app.locs]=findpeaks(abs(app.stackMFN),'MINPEAKHEIGHT',app.MinPeak); 
            app.GainFunc = ones(length(app.stackMF),1);
            for i = 1:length(app.stackMF)
                if app.t(i) < app.TimeDelay
                    app.GainFunc(i) = 1;
                elseif app.t(i) >= app.TimeDelay  && app.t(app.EndSample) >= app.t(i)
                    app.GainFunc(i) = exp(app.K*(app.t(i) - app.TimeDelay));
                else app.t(i) >= app.t(app.EndSample)
                    app.GainFunc(i) = app.ExpGain*exp(-app.K*(app.t(i) - app.EndTime));
                end
            end
            app.stackMFG = (app.GainFunc.*app.stackMF);
            app.stackMFGN = app.stackMFG/max(abs(app.stackMFG));
            [app.pks,app.locs]=findpeaks(abs(app.stackMFGN),'MINPEAKHEIGHT',app.MinPeak);
            OnOff = app.FilterSwitch.Value;


                switch OnOff
                    case 'Yes'
                        
                        delete(plot(app.R3, app.x, app.stack, '-b'));
                        delete(plot(app.D3, app.fB,10*log10(app.pxxB),'-', Color='b'));   
                        [app.pxxB,app.fB,app.pxxcB] = periodogram(app.stackMFGN,rectwin(length(app.stack)),length(app.stackMFGN),app.fs);
                        plot(app.R3, app.x,app.stackMFGN, '-b', app.x(app.locs([1,end])), app.stackMFGN(app.locs([1,end])), 'sr');
                        title(app.R3, 'Stacked and filtered record reflectogram');
                        xlim(app.R3,[app.XMin, app.XMax]);

                        app.N = num2cell(round(app.x(app.locs([1,end]))*100)./100);
                        text(app.R3, app.x(app.locs([1,end]))+app.offsetX, app.stackMFN(app.locs([1, end]))-app.offsetY, app.N , ...
                        'Color','red','FontName', 'Candara', 'FontSize',12);
                        [app.pxxC,app.fC,app.pxxcC] = periodogram(app.stackMFGN,rectwin(length(app.stackMFGN)),length(app.stackMFGN),app.fs);
                        plot(app.D3, app.fB,10*log10(app.pxxC),'-', Color='b');
                        title(app.D3, 'Spectral density of the stacked and filtered record');
                                    
                    case 'No'
                        delete( plot(app.R3, app.x,app.stackMFGN, '-b', app.x(app.locs([1,end])), app.stackMFGN(app.locs([1,end])), 'sr'));
                        delete( plot(app.D3, app.fB,10*log10(app.pxxC),'-', Color='b'));
                        plot(app.R3, app.x, app.stack, '-b');   %Grafica individual  
                        title (app.R3, 'Stacked Record Reflectogram');
                        [app.pxxB,app.fB,app.pxxcB] = periodogram(app.stack,rectwin(length(app.stack)),length(app.stack),app.fs);
                        plot(app.D3, app.fB,10*log10(app.pxxB),'-', Color='b');
                        title (app.D3, 'Spectral Density of the Stacked Record');

                end
        end

        % Button pushed function: PrintReportButton
        function PrintReportButtonPushed(app, event)
            app.NumFile = app.RecordNSpinner.Value;
            app.DataM = app.data(:,app.NumFile);
            app.date = app.DateDatePicker.Value;
            app.DateString = datestr(app.date );
            app.stackM = app.stack/size(app.data,2);
            app.stackM = app.stackM/max(abs([min(app.stackM) max(app.stackM)]));
            app.PileName = app.PileNameEditField.Value;
            app.velocity = app.VcmsEditField.Value;                 %Concrete velocity m/s
            app.fs = app.fsHzEditField.Value;                       %frequency sampling in Hz
            app.PileLength = app.LmEditField.Value;                 %in meters
            app.Delay = app.DelayEditField.Value;
            app.ExpGain = app.ExpGainEditField.Value;
            app.Order = 1;
            app.LCF = app.LCFHzEditField.Value;
            app.HCF = app.HCFHzEditField.Value;
            app.date = app.DateDatePicker.Value;
            app.MinPeak = 0.75;
            app.offsetX = 0.1;
            app.offsetY = 0.1;
            app.XMin = 0.0;
            app.XMax = 1.2*app.PileLength;
            app.YMin = -1.0;
            app.YMax = 1.0;
            app.dt = 1/app.fs;
            app.dx = app.dt*app.velocity/2;
            app.t = 0:app.dt:(length(app.stack)-1)*app.dt;
            app.x = 0:app.dx:(length(app.stack)-1)*app.dx;
            [app.b,app.a] = butter(app.Order,[app.LCF app.HCF]/(app.fs/2),'bandpass');
            app.stackMF = filter(app.b,app.a,app.stackM);
            app.EndLength = 1.10*app.PileLength; %in meters
            app.EndSample = floor((app.EndLength)/app.dx);
            app.EndTime = app.EndSample*app.dt;
            app.TimeDelay = (app.Delay/100)*app.dt*app.EndSample;
            app.XDelay = (app.Delay/100)*app.dx*app.EndSample;
            app.K = log(app.ExpGain)/((app.EndSample*app.dt)-app.TimeDelay);
            app.GainFunc = ones(length(app.stackMF),1);
            for i = 1:length(app.stackMF)
                if app.t(i) < app.TimeDelay
                    app.GainFunc(i) = 1;
                elseif app.t(i) >= app.TimeDelay  && app.t(app.EndSample) >= app.t(i)
                    app.GainFunc(i) = exp(app.K*(app.t(i) - app.TimeDelay));
                else app.t(i) >= app.t(app.EndSample)
                    app.GainFunc(i) = app.ExpGain*exp(-app.K*(app.t(i) - app.EndTime));
                end
            end
            app.stackMFN = app.stackMF/max(abs(app.stackMF));
            [app.pks,app.locs]=findpeaks(abs(app.stackMFN),'MINPEAKHEIGHT',app.MinPeak); 
            app.stackMFG = (app.GainFunc.*app.stackMF);
            app.stackMFGN = app.stackMFG/max(abs(app.stackMFG));
            [app.pks,app.locs]=findpeaks(abs(app.stackMFGN),'MINPEAKHEIGHT',app.MinPeak);
            app.figout = figure('visible','off');
            plot(app.x,app.stackMFGN, '-b', app.x(app.locs([1,end])), app.stackMFGN(app.locs([1,end])), 'sr');
            app.N = num2cell(round(app.x(app.locs([1,end]))*100)./100);
            text(app.x(app.locs([1,end]))+app.offsetX, app.stackMFGN(app.locs([1, end]))-app.offsetY, app.N , ...
                'Color','red','FontName', 'Helvetica', 'FontSize',8);
            text(app.XMin,app.YMin-0.4, strcat('Fecha: ',app.DateString),'Color',[0.3 0.3 0.3],'FontName', ...
                'Helvetica', 'FontSize',7);
            text(app.XMin,app.YMin-0.55, sprintf('Pilote: %s', app.PileName),'Color',[0.3 0.3 0.3],'FontName', ...
                'Helvetica', 'FontSize',7);
            text(app.XMin,app.YMin-0.7, sprintf('G = %d', app.ExpGain),'Color',[0.3 0.3 0.3], ...
                'FontName','Helvetica', 'FontSize',7);
            text(app.XMin,app.YMin-0.85, sprintf('Retardo = %2.2f m', app.XDelay),'Color',[0.3 0.3 0.3], ...
                'FontName','Helvetica', 'FontSize',7);
            text(app.XMin,app.YMin-1.0, sprintf('Long. Est. = %2.2f m', app.x(app.locs(end)) - app.x(app.locs(1))), ...
                'Color',[0.3 0.3 0.3],'FontName','Helvetica', 'FontSize',7);
            text(app.XMin,app.YMin-1.15, sprintf('Procesado en Easy PIT Pro - Subsuelo 3D SAS, 2022'), ...
                'Color',[0 0 0],'FontName','Helvetica', 'FontSize',7);
            set(gcf,'color','w');
            title('Reflectograma PIT', 'FontName','Helvetica','FontSize', 9, 'Color', [0.3 0.3 0.3]);
            xlabel('Longitud (m)', 'FontName','Helvetica','FontSize', 8);
            ylabel('Amplitud normalizada', 'FontName','Helvetica','FontSize', 8);
            set(gca, 'YGrid', 'on', 'XGrid', 'on', 'XMinorGrid', 'on', 'YMinorGrid', 'on', ...
                'XColor', [0.3 0.3 0.3], 'YColor',[0.3 0.3 0.3], ...
                'FontName','Helvetica','FontSize',8, 'TickDir', 'out');
            xlim([app.XMin, app.XMax]);
            ylim([app.YMin, app.YMax]);
            pbaspect([1 1/3 1]);

            app.tipo = {
                '*.jpg';'*.png';'*.tif';'*.pdf';'*.eps'};
            [filename,filepath] = uiputfile(app.tipo);
            if ischar(filename)
                saveas(app.figout,[filepath filename]);
            end
%             print(app.figout, '-dpdf', strcat(app.PileName,'.pdf'), '-r3000')
        end

        % Value changed function: LmEditField
        function LmEditFieldValueChanged(app, event)
            app.samplenum = length(app.stack);
            app.recordnum = size(app.data,2);
            app.stackM = app.stack/size(app.data,2);
            app.stackM = app.stackM/max(abs([min(app.stackM) max(app.stackM)]));
            app.PileName = app.PileNameEditField.Value;
            app.velocity = app.VcmsEditField.Value;                 %Concrete velocity m/s
            app.fs = app.fsHzEditField.Value;                       %frequency sampling in Hz
            app.PileLength = app.LmEditField.Value;                 %in meters
            app.Delay = app.DelayEditField.Value;
            app.ExpGain = app.ExpGainEditField.Value;
            app.Order = 1;
            app.LCF = app.LCFHzEditField.Value;
            app.HCF = app.HCFHzEditField.Value;
            app.date = app.DateDatePicker.Value;
            app.MinPeak = 0.75;
            app.offsetX = 0.1;
            app.offsetY = 0.1;
            app.XMin = 0.0;
            app.XMax = 1.2*app.PileLength;
            app.YMin = -1.0;
            app.YMax = 1.0;
            app.dt = 1/app.fs;
            app.dx = app.dt*app.velocity/2;
            app.t = 0:app.dt:(length(app.stack)-1)*app.dt;
            app.x = 0:app.dx:(length(app.stack)-1)*app.dx;
            [app.b,app.a] = butter(app.Order,[app.LCF app.HCF]/(app.fs/2),'bandpass');
            app.stackMF = filter(app.b,app.a,app.stackM);
            app.EndLength = 1.10*app.PileLength; %in meters
            app.EndSample = floor((app.EndLength)/app.dx);
            app.EndTime = app.EndSample*app.dt;
            app.TimeDelay = (app.Delay/100)*app.dt*app.EndSample;
            app.XDelay = (app.Delay/100)*app.dx*app.EndSample;
            app.K = log(app.ExpGain)/((app.EndSample*app.dt)-app.TimeDelay);
            app.GainFunc = ones(length(app.stackMF),1);
            for i = 1:length(app.stackMF)
                if app.t(i) < app.TimeDelay
                    app.GainFunc(i) = 1;
                elseif app.t(i) >= app.TimeDelay  && app.t(app.EndSample) >= app.t(i)
                    app.GainFunc(i) = exp(app.K*(app.t(i) - app.TimeDelay));
                else app.t(i) >= app.t(app.EndSample)
                    app.GainFunc(i) = app.ExpGain*exp(-app.K*(app.t(i) - app.EndTime));
                end
            end
            app.stackMFG = (app.GainFunc.*app.stackMF);
            app.stackMFGN = app.stackMFG/max(abs(app.stackMFG));
            [app.pks,app.locs]=findpeaks(abs(app.stackMFGN),'MINPEAKHEIGHT',app.MinPeak); %Graficas
            %----Individuales----
            plot(app.R1, app.x,app.DataM, '-b');
            title(app.R1, sprintf('Individual Record Reflectogram. N°= %d', app.NumFile));
            xlim(app.R1,[app.XMin, app.XMax]);
            [app.pxx,app.f,app.pxxc] = periodogram(app.DataM,rectwin(length(app.DataM)),length(app.DataM),app.fs,...
            'ConfidenceLevel',0.95);
            plot(app.D1, app.f,10*log10(app.pxx),'-', Color='b');
             title(app.D1, sprintf('Individual Record Spectral Density. N°= %d', app.NumFile));
            %----Ganancias----
            plot(app.Gain, app.x, app.GainFunc, '-b');

            %----Apilado---
            plot(app.R3, app.x, app.stack, '-b');   %Grafica individual       
            [app.pxxB,app.fB,app.pxxcB] = periodogram(app.stack,rectwin(length(app.stack)),length(app.stack),app.fs);
            plot(app.D3, app.fB,10*log10(app.pxxB),'-', Color='b');
           
        end

        % Value changed function: LCFHzEditField
        function LCFHzEditFieldValueChanged(app, event)
            app.NumFile = app.RecordNSpinner.Value;
            app.DataM = app.data(:,app.NumFile);
            app.date = app.DateDatePicker.Value;
            app.DateString = datestr(app.date );
            app.stackM = app.stack/size(app.data,2);
            app.stackM = app.stackM/max(abs([min(app.stackM) max(app.stackM)]));
            app.PileName = app.PileNameEditField.Value;
            app.velocity = app.VcmsEditField.Value;                 %Concrete velocity m/s
            app.fs = app.fsHzEditField.Value;                       %frequency sampling in Hz
            app.PileLength = app.LmEditField.Value;                 %in meters
            app.Delay = app.DelayEditField.Value;
            app.ExpGain = app.ExpGainEditField.Value;
            app.Order = 1;
            app.LCF = app.LCFHzEditField.Value;
            app.HCF = app.HCFHzEditField.Value;
            app.date = app.DateDatePicker.Value;
            app.MinPeak = 0.75;
            app.offsetX = 0.1;
            app.offsetY = 0.1;
            app.XMin = 0.0;
            app.XMax = 1.2*app.PileLength;
            app.YMin = -1.0;
            app.YMax = 1.0;
            app.dt = 1/app.fs;
            app.dx = app.dt*app.velocity/2;
            app.t = 0:app.dt:(length(app.stack)-1)*app.dt;
            app.x = 0:app.dx:(length(app.stack)-1)*app.dx;
            [app.b,app.a] = butter(app.Order,[app.LCF app.HCF]/(app.fs/2),'bandpass');
            app.stackMF = filter(app.b,app.a,app.stackM);
            app.EndLength = 1.10*app.PileLength; %in meters
            app.EndSample = floor((app.EndLength)/app.dx);
            app.EndTime = app.EndSample*app.dt;
            app.TimeDelay = (app.Delay/100)*app.dt*app.EndSample;
            app.XDelay = (app.Delay/100)*app.dx*app.EndSample;
            app.K = log(app.ExpGain)/((app.EndSample*app.dt)-app.TimeDelay);
            app.GainFunc = ones(length(app.stackMF),1);
            for i = 1:length(app.stackMF)
                if app.t(i) < app.TimeDelay
                    app.GainFunc(i) = 1;
                elseif app.t(i) >= app.TimeDelay  && app.t(app.EndSample) >= app.t(i)
                    app.GainFunc(i) = exp(app.K*(app.t(i) - app.TimeDelay));
                else app.t(i) >= app.t(app.EndSample)
                    app.GainFunc(i) = app.ExpGain*exp(-app.K*(app.t(i) - app.EndTime));
                end
            end
            app.stackMFN = app.stackMF/max(abs(app.stackMF));
            [app.pks,app.locs]=findpeaks(abs(app.stackMFN),'MINPEAKHEIGHT',app.MinPeak); 
            %Graficas
            [app.pxxB,app.fB,app.pxxcB] = periodogram(app.stackMFGN,rectwin(length(app.stack)),length(app.stackMFGN),app.fs);
            plot(app.R3, app.x,app.stackMFGN, '-b', app.x(app.locs([1,end])), app.stackMFGN(app.locs([1,end])), 'sr');
            title(app.R3, 'Stacked and filtered record reflectogram');

            app.N = num2cell(round(app.x(app.locs([1,end]))*100)./100);
            text(app.R3, app.x(app.locs([1,end]))+app.offsetX, app.stackMFN(app.locs([1, end]))-app.offsetY, app.N , ...
            'Color','red','FontName', 'Candara', 'FontSize',12);
            [app.pxxC,app.fC,app.pxxcC] = periodogram(app.stackMFGN,rectwin(length(app.stackMFGN)),length(app.stackMFGN),app.fs);
            plot(app.D3, app.fB,10*log10(app.pxxC),'-', Color='b');
            title(app.D3, 'Spectral density of the stacked and filtered record');

            %----Ganancias----
            plot(app.Gain, app.x, app.GainFunc, '-b');

        end

        % Value changed function: HCFHzEditField
        function HCFHzEditFieldValueChanged(app, event)
            app.samplenum = length(app.stack);
            app.recordnum = size(app.data,2);
            app.stackM = app.stack/size(app.data,2);
            app.stackM = app.stackM/max(abs([min(app.stackM) max(app.stackM)]));
            app.PileName = app.PileNameEditField.Value;
            app.velocity = app.VcmsEditField.Value;                 %Concrete velocity m/s
            app.fs = app.fsHzEditField.Value;                       %frequency sampling in Hz
            app.PileLength = app.LmEditField.Value;                 %in meters
            app.Delay = app.DelayEditField.Value;
            app.ExpGain = app.ExpGainEditField.Value;
            app.Order = 1;
            app.LCF = app.LCFHzEditField.Value;
            app.HCF = app.HCFHzEditField.Value;
            app.date = app.DateDatePicker.Value;
            app.MinPeak = 0.75;
            app.offsetX = 0.1;
            app.offsetY = 0.1;
            app.XMin = 0.0;
            app.XMax = 1.2*app.PileLength;
            app.YMin = -1.0;
            app.YMax = 1.0;
            app.dt = 1/app.fs;
            app.dx = app.dt*app.velocity/2;
            app.t = 0:app.dt:(length(app.stack)-1)*app.dt;
            app.x = 0:app.dx:(length(app.stack)-1)*app.dx;
            [app.b,app.a] = butter(app.Order,[app.LCF app.HCF]/(app.fs/2),'bandpass');
            app.stackMF = filter(app.b,app.a,app.stackM);
            app.EndLength = 1.10*app.PileLength; %in meters
            app.EndSample = floor((app.EndLength)/app.dx);
            app.EndTime = app.EndSample*app.dt;
            app.TimeDelay = (app.Delay/100)*app.dt*app.EndSample;
            app.XDelay = (app.Delay/100)*app.dx*app.EndSample;
            app.K = log(app.ExpGain)/((app.EndSample*app.dt)-app.TimeDelay);
            app.GainFunc = ones(length(app.stackMF),1);
            for i = 1:length(app.stackMF)
                if app.t(i) < app.TimeDelay
                    app.GainFunc(i) = 1;
                elseif app.t(i) >= app.TimeDelay  && app.t(app.EndSample) >= app.t(i)
                    app.GainFunc(i) = exp(app.K*(app.t(i) - app.TimeDelay));
                else app.t(i) >= app.t(app.EndSample)
                    app.GainFunc(i) = app.ExpGain*exp(-app.K*(app.t(i) - app.EndTime));
                end
            end
            app.stackMFG = (app.GainFunc.*app.stackMF);
            app.stackMFGN = app.stackMFG/max(abs(app.stackMFG));
            [app.pks,app.locs]=findpeaks(abs(app.stackMFGN),'MINPEAKHEIGHT',app.MinPeak); %Graficas
            %Graficas
            [app.pxxB,app.fB,app.pxxcB] = periodogram(app.stackMFGN,rectwin(length(app.stack)),length(app.stackMFGN),app.fs);
            plot(app.R3, app.x,app.stackMFGN, '-b', app.x(app.locs([1,end])), app.stackMFGN(app.locs([1,end])), 'sr');
            title(app.R3, 'Stacked and filtered record reflectogram');
        
            app.N = num2cell(round(app.x(app.locs([1,end]))*100)./100);
            text(app.R3, app.x(app.locs([1,end]))+app.offsetX, app.stackMFN(app.locs([1, end]))-app.offsetY, app.N , ...
            'Color','red','FontName', 'Candara', 'FontSize',12);
            [app.pxxC,app.fC,app.pxxcC] = periodogram(app.stackMFGN,rectwin(length(app.stackMFGN)),length(app.stackMFGN),app.fs);
            plot(app.D3, app.fB,10*log10(app.pxxC),'-', Color='b');
            title(app.D3, 'Spectral density of the stacked and filtered record');

            
        end

        % Value changed function: DelayEditField
        function DelayEditFieldValueChanged(app, event)
            app.samplenum = length(app.stack);
            app.recordnum = size(app.data,2);
            app.stackM = app.stack/size(app.data,2);
            app.stackM = app.stackM/max(abs([min(app.stackM) max(app.stackM)]));
            app.PileName = app.PileNameEditField.Value;
            app.velocity = app.VcmsEditField.Value;                 %Concrete velocity m/s
            app.fs = app.fsHzEditField.Value;                       %frequency sampling in Hz
            app.PileLength = app.LmEditField.Value;                 %in meters
            app.Delay = app.DelayEditField.Value;
            app.ExpGain = app.ExpGainEditField.Value;
            app.Order = 1;
            app.LCF = app.LCFHzEditField.Value;
            app.HCF = app.HCFHzEditField.Value;
            app.date = app.DateDatePicker.Value;
            app.MinPeak = 0.75;
            app.offsetX = 0.1;
            app.offsetY = 0.1;
            app.XMin = 0.0;
            app.XMax = 1.2*app.PileLength;
            app.YMin = -1.0;
            app.YMax = 1.0;
            app.dt = 1/app.fs;
            app.dx = app.dt*app.velocity/2;
            app.t = 0:app.dt:(length(app.stack)-1)*app.dt;
            app.x = 0:app.dx:(length(app.stack)-1)*app.dx;
            [app.b,app.a] = butter(app.Order,[app.LCF app.HCF]/(app.fs/2),'bandpass');
            app.stackMF = filter(app.b,app.a,app.stackM);
            app.EndLength = 1.10*app.PileLength; %in meters
            app.EndSample = floor((app.EndLength)/app.dx);
            app.EndTime = app.EndSample*app.dt;
            app.TimeDelay = (app.Delay/100)*app.dt*app.EndSample;
            app.XDelay = (app.Delay/100)*app.dx*app.EndSample;
            app.K = log(app.ExpGain)/((app.EndSample*app.dt)-app.TimeDelay);
            app.GainFunc = ones(length(app.stackMF),1);
            for i = 1:length(app.stackMF)
                if app.t(i) < app.TimeDelay
                    app.GainFunc(i) = 1;
                elseif app.t(i) >= app.TimeDelay  && app.t(app.EndSample) >= app.t(i)
                    app.GainFunc(i) = exp(app.K*(app.t(i) - app.TimeDelay));
                else app.t(i) >= app.t(app.EndSample)
                    app.GainFunc(i) = app.ExpGain*exp(-app.K*(app.t(i) - app.EndTime));
                end
            end
            app.stackMFG = (app.GainFunc.*app.stackMF);
            app.stackMFGN = app.stackMFG/max(abs(app.stackMFG));
            [app.pks,app.locs]=findpeaks(abs(app.stackMFGN),'MINPEAKHEIGHT',app.MinPeak); %Graficas
            %Graficas 
            [app.pxxB,app.fB,app.pxxcB] = periodogram(app.stackMFGN,rectwin(length(app.stack)),length(app.stackMFGN),app.fs);
            plot(app.R3, app.x,app.stackMFGN, '-b', app.x(app.locs([1,end])), app.stackMFGN(app.locs([1,end])), 'sr');
            title(app.R3, 'Stacked and filtered record reflectogram');

            app.N = num2cell(round(app.x(app.locs([1,end]))*100)./100);
            text(app.R3, app.x(app.locs([1,end]))+app.offsetX, app.stackMFN(app.locs([1, end]))-app.offsetY, app.N , ...
            'Color','red','FontName', 'Candara', 'FontSize',12);
            [app.pxxC,app.fC,app.pxxcC] = periodogram(app.stackMFGN,rectwin(length(app.stackMFGN)),length(app.stackMFGN),app.fs);
            plot(app.D3, app.fB,10*log10(app.pxxC),'-', Color='b');
            title(app.D3, 'Spectral density of the stacked and filtered record');
            %----Ganancias----
            plot(app.Gain, app.x, app.GainFunc, '-b');

        end

        % Value changed function: ExpGainEditField
        function ExpGainEditFieldValueChanged(app, event)
            app.samplenum = length(app.stack);
            app.recordnum = size(app.data,2);
            app.stackM = app.stack/size(app.data,2);
            app.stackM = app.stackM/max(abs([min(app.stackM) max(app.stackM)]));
            app.PileName = app.PileNameEditField.Value;
            app.velocity = app.VcmsEditField.Value;                 %Concrete velocity m/s
            app.fs = app.fsHzEditField.Value;                       %frequency sampling in Hz
            app.PileLength = app.LmEditField.Value;                 %in meters
            app.Delay = app.DelayEditField.Value;
            app.ExpGain = app.ExpGainEditField.Value;
            app.Order = 1;
            app.LCF = app.LCFHzEditField.Value;
            app.HCF = app.HCFHzEditField.Value;
            app.date = app.DateDatePicker.Value;
            app.MinPeak = 0.75;
            app.offsetX = 0.1;
            app.offsetY = 0.1;
            app.XMin = 0.0;
            app.XMax = 1.2*app.PileLength;
            app.YMin = -1.0;
            app.YMax = 1.0;
            app.dt = 1/app.fs;
            app.dx = app.dt*app.velocity/2;
            app.t = 0:app.dt:(length(app.stack)-1)*app.dt;
            app.x = 0:app.dx:(length(app.stack)-1)*app.dx;
            [app.b,app.a] = butter(app.Order,[app.LCF app.HCF]/(app.fs/2),'bandpass');
            app.stackMF = filter(app.b,app.a,app.stackM);
            app.EndLength = 1.10*app.PileLength; %in meters
            app.EndSample = floor((app.EndLength)/app.dx);
            app.EndTime = app.EndSample*app.dt;
            app.TimeDelay = (app.Delay/100)*app.dt*app.EndSample;
            app.XDelay = (app.Delay/100)*app.dx*app.EndSample;
            app.K = log(app.ExpGain)/((app.EndSample*app.dt)-app.TimeDelay);
            app.GainFunc = ones(length(app.stackMF),1);
            for i = 1:length(app.stackMF)
                if app.t(i) < app.TimeDelay
                    app.GainFunc(i) = 1;
                elseif app.t(i) >= app.TimeDelay  && app.t(app.EndSample) >= app.t(i)
                    app.GainFunc(i) = exp(app.K*(app.t(i) - app.TimeDelay));
                else app.t(i) >= app.t(app.EndSample)
                    app.GainFunc(i) = app.ExpGain*exp(-app.K*(app.t(i) - app.EndTime));
                end
            end
            app.stackMFG = (app.GainFunc.*app.stackMF);
            app.stackMFGN = app.stackMFG/max(abs(app.stackMFG));
            [app.pks,app.locs]=findpeaks(abs(app.stackMFGN),'MINPEAKHEIGHT',app.MinPeak); %Graficas
            [app.pxxB,app.fB,app.pxxcB] = periodogram(app.stackMFGN,rectwin(length(app.stack)),length(app.stackMFGN),app.fs);
            plot(app.R3, app.x,app.stackMFGN, '-b', app.x(app.locs([1,end])), app.stackMFGN(app.locs([1,end])), 'sr');
            title(app.R3, 'Stacked and filtered record reflectogram');

            app.N = num2cell(round(app.x(app.locs([1,end]))*100)./100);
            text(app.R3, app.x(app.locs([1,end]))+app.offsetX, app.stackMFN(app.locs([1, end]))-app.offsetY, app.N , ...
            'Color','red','FontName', 'Candara', 'FontSize',12);
            [app.pxxC,app.fC,app.pxxcC] = periodogram(app.stackMFGN,rectwin(length(app.stackMFGN)),length(app.stackMFGN),app.fs);
            plot(app.D3, app.fB,10*log10(app.pxxC),'-', Color='b');
            title(app.D3, 'Spectral density of the stacked and filtered record');
            %----Ganancias----
            plot(app.Gain, app.x, app.GainFunc, '-b');


            app.o = transpose(app.FileName);
            app.o = array2table(app.o );
            app.UITable.Data = app.o;
           

        end

        % Button pushed function: StackButton
        function StackButtonPushed(app, event)
            app.FigoutStack = figure('visible','off'); 
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 280 363];
            app.UIFigure.Name = 'MATLAB App';

            % Create UITable
            app.UITable = uitable(app.UIFigure);
            app.UITable.BackgroundColor = [0.9412 0.9412 0.9412;0.8 0.8 0.8];
            app.UITable.ColumnName = {'File Name'};
            app.UITable.RowName = {'1'; '2'; '3'; '4'; '5'; '6'; '7'; '8'; '9'; '10'; '11'; '12'; ''};
            app.UITable.ColumnEditable = true;
            app.UITable.RowStriping = 'off';
            app.UITable.FontName = 'Candara';
            app.UITable.Position = [20 51 244 270];

            % Create OKButton
            app.OKButton = uibutton(app.UIFigure, 'push');
            app.OKButton.ButtonPushedFcn = createCallbackFcn(app, @OKButtonPushed, true);
            app.OKButton.FontName = 'Candara';
            app.OKButton.Position = [91 22 100 22];
            app.OKButton.Text = 'OK';

            % Create Step4Label
            app.Step4Label = uilabel(app.UIFigure);
            app.Step4Label.BackgroundColor = [0.9412 0.9412 0.9412];
            app.Step4Label.FontName = 'Candara';
            app.Step4Label.FontSize = 14;
            app.Step4Label.Tooltip = {'Choose records'};
            app.Step4Label.Position = [18 334 43 22];
            app.Step4Label.Text = 'Step 4';

            app.o = transpose(app.FileName);
            app.o = cell2table(app.o );
            app.UITable.Data = app.o;
            app.UIFigure.Visible = 'on';

            function OKButtonPushed(app, event)
                [m,n] = size (app.UITable.Data);
                app.UITableData = app.UITable.Data;
                app.d = cell(m+1,n);
                app.d(1:m,1:n)= table2cell(app.UITableData);
                app.FileName2 = app.d(1:m,1:n);
                delete(app.UIFigure)
                
                %# find empty cells
                app.emptyCells = cellfun(@isempty,app.FileName2);
                %# remove empty cells
                app.FileName2(app.emptyCells) = [];
                
                if length(app.d(1:m,1:n)) < 2 || length(strfind(app.FileName2, '.TXT')) < 2 || length(strfind(app.FileName2, '.txt')) < 2
                app.data = load(app.FileName2);
                app.stack = app.data;
                    if size(app.data,2) > size(app.data,1)
                    app.stack = app.stack';
                    end
                    app.NumFiles = 1;
                    app.k = 1;
                else
                    app.k = 1; 
                    app.NumFiles = length(app.FileName2);
                    app.data = zeros(length(load(app.FileName2{app.k})), length(app.FileName2));
                    app.data(:,app.k) = load(app.FileName2{app.k});
                    app.stack = app.data(:,app.k);
                    for k = 2:length(app.FileName2)
                        app.data(:,k) = load(app.FileName2{k});
                        app.stack = app.stack + app.data(:,k);
                        waitbar(k/length(app.FileName2));
                    end
                end
                plot(app.R3, app.x, app.stack, '-b');   %Grafica individual  
                title (app.R3, 'Stacked Record Reflectogram');
                [app.pxxB,app.fB,app.pxxcB] = periodogram(app.stack,rectwin(length(app.stack)),length(app.stack),app.fs);
                plot(app.D3, app.fB,10*log10(app.pxxB),'-', Color='b');
                title (app.D3, 'Spectral Density of the Stacked Record');
            end  
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Color = [0.9412 0.9412 0.9412];
            app.UIFigure.Position = [100 50 1143 660];
            app.UIFigure.Name = 'MATLAB App';

            % Create Gain
            app.Gain = uiaxes(app.UIFigure);
            title(app.Gain, 'Exponential Gain Function')
            xlabel(app.Gain, 'Length (m)')
            ylabel(app.Gain, 'Gain')
            zlabel(app.Gain, 'Z')
            app.Gain.FontName = 'Candara';
            app.Gain.XGrid = 'on';
            app.Gain.XMinorGrid = 'on';
            app.Gain.YGrid = 'on';
            app.Gain.YMinorGrid = 'on';
            app.Gain.FontSize = 14;
            app.Gain.Position = [460 262 400 178];

            % Create D1
            app.D1 = uiaxes(app.UIFigure);
            title(app.D1, 'Individual Record Spectral Density')
            xlabel(app.D1, 'Frequency (Hz)')
            ylabel(app.D1, 'Spectral Density (dB/Hz)')
            zlabel(app.D1, 'Z')
            app.D1.FontName = 'Candara';
            app.D1.XMinorGrid = 'on';
            app.D1.YGrid = 'on';
            app.D1.YMinorGrid = 'on';
            app.D1.FontSize = 14;
            app.D1.Position = [680 450 400 178];

            % Create R1
            app.R1 = uiaxes(app.UIFigure);
            title(app.R1, 'Individual Record Reflectogram')
            xlabel(app.R1, 'Length (m)')
            ylabel(app.R1, 'Amplitude')
            zlabel(app.R1, 'Z')
            app.R1.FontName = 'Candara';
            app.R1.XGrid = 'on';
            app.R1.XMinorGrid = 'on';
            app.R1.YGrid = 'on';
            app.R1.YMinorGrid = 'on';
            app.R1.FontSize = 14;
            app.R1.MinorGridAlpha = 0.25;
            app.R1.Position = [237 450 400 178];

            % Create R3
            app.R3 = uiaxes(app.UIFigure);
            title(app.R3, 'Stacked Record Reflectogram')
            xlabel(app.R3, 'Length (m)')
            ylabel(app.R3, 'Normalized Amplitude')
            zlabel(app.R3, 'Z')
            app.R3.FontName = 'Candara';
            app.R3.XGrid = 'on';
            app.R3.XMinorGrid = 'on';
            app.R3.YGrid = 'on';
            app.R3.YMinorGrid = 'on';
            app.R3.FontSize = 14;
            app.R3.Position = [237 77 400 178];

            % Create D3
            app.D3 = uiaxes(app.UIFigure);
            title(app.D3, 'Spectral Density of the Stacked Record')
            xlabel(app.D3, 'Frequency (Hz)')
            ylabel(app.D3, 'Spectral Density (dB/Hz)')
            zlabel(app.D3, 'Z')
            app.D3.FontName = 'Candara';
            app.D3.XMinorGrid = 'on';
            app.D3.YGrid = 'on';
            app.D3.YMinorGrid = 'on';
            app.D3.FontSize = 14;
            app.D3.Position = [680 77 400 178];

            % Create Panel
            app.Panel = uipanel(app.UIFigure);
            app.Panel.Tooltip = {'Choose records'};
            app.Panel.FontName = 'Candara';
            app.Panel.FontWeight = 'bold';
            app.Panel.FontSize = 14;
            app.Panel.Position = [25 16 198 632];

            % Create PileNameEditFieldLabel
            app.PileNameEditFieldLabel = uilabel(app.Panel);
            app.PileNameEditFieldLabel.HorizontalAlignment = 'right';
            app.PileNameEditFieldLabel.FontName = 'Candara';
            app.PileNameEditFieldLabel.FontSize = 14;
            app.PileNameEditFieldLabel.Tooltip = {'Enter name of the pile and parameters of speed, frequency and length'};
            app.PileNameEditFieldLabel.Position = [22 514 65 22];
            app.PileNameEditFieldLabel.Text = 'Pile Name';

            % Create PileNameEditField
            app.PileNameEditField = uieditfield(app.Panel, 'text');
            app.PileNameEditField.HorizontalAlignment = 'right';
            app.PileNameEditField.FontName = 'Candara';
            app.PileNameEditField.FontSize = 14;
            app.PileNameEditField.Position = [98 514 67 22];
            app.PileNameEditField.Value = 'Pile 0';

            % Create VcmsEditFieldLabel
            app.VcmsEditFieldLabel = uilabel(app.Panel);
            app.VcmsEditFieldLabel.HorizontalAlignment = 'right';
            app.VcmsEditFieldLabel.FontName = 'Candara';
            app.VcmsEditFieldLabel.FontSize = 14;
            app.VcmsEditFieldLabel.Tooltip = {'Enter name of the pile and parameters of speed, frequency and length'};
            app.VcmsEditFieldLabel.Position = [34 484 50 22];
            app.VcmsEditFieldLabel.Text = 'Vc(m/s)';

            % Create VcmsEditField
            app.VcmsEditField = uieditfield(app.Panel, 'numeric');
            app.VcmsEditField.ValueChangedFcn = createCallbackFcn(app, @VcmsEditFieldValueChanged, true);
            app.VcmsEditField.FontName = 'Candara';
            app.VcmsEditField.FontSize = 14;
            app.VcmsEditField.Position = [98 484 67 22];
            app.VcmsEditField.Value = 3352;

            % Create fsHzEditFieldLabel
            app.fsHzEditFieldLabel = uilabel(app.Panel);
            app.fsHzEditFieldLabel.HorizontalAlignment = 'right';
            app.fsHzEditFieldLabel.FontName = 'Candara';
            app.fsHzEditFieldLabel.FontSize = 14;
            app.fsHzEditFieldLabel.Tooltip = {'Enter name of the pile and parameters of speed, frequency and length'};
            app.fsHzEditFieldLabel.Position = [43 454 41 22];
            app.fsHzEditFieldLabel.Text = 'fs(Hz)';

            % Create fsHzEditField
            app.fsHzEditField = uieditfield(app.Panel, 'numeric');
            app.fsHzEditField.ValueDisplayFormat = '%.0f';
            app.fsHzEditField.ValueChangedFcn = createCallbackFcn(app, @fsHzEditFieldValueChanged, true);
            app.fsHzEditField.FontName = 'Candara';
            app.fsHzEditField.FontSize = 14;
            app.fsHzEditField.BackgroundColor = [1 1 0];
            app.fsHzEditField.Position = [98 454 67 22];
            app.fsHzEditField.Value = 52734;

            % Create LmEditFieldLabel
            app.LmEditFieldLabel = uilabel(app.Panel);
            app.LmEditFieldLabel.HorizontalAlignment = 'right';
            app.LmEditFieldLabel.FontName = 'Candara';
            app.LmEditFieldLabel.FontSize = 14;
            app.LmEditFieldLabel.Position = [50 424 34 22];
            app.LmEditFieldLabel.Text = 'L(m)';

            % Create LmEditField
            app.LmEditField = uieditfield(app.Panel, 'numeric');
            app.LmEditField.ValueChangedFcn = createCallbackFcn(app, @LmEditFieldValueChanged, true);
            app.LmEditField.FontName = 'Candara';
            app.LmEditField.FontSize = 14;
            app.LmEditField.Position = [98 424 67 22];
            app.LmEditField.Value = 20;

            % Create RecordNSpinnerLabel
            app.RecordNSpinnerLabel = uilabel(app.Panel);
            app.RecordNSpinnerLabel.HorizontalAlignment = 'right';
            app.RecordNSpinnerLabel.FontName = 'Candara';
            app.RecordNSpinnerLabel.FontSize = 14;
            app.RecordNSpinnerLabel.Position = [21 364 64 22];
            app.RecordNSpinnerLabel.Text = 'Record N°';

            % Create RecordNSpinner
            app.RecordNSpinner = uispinner(app.Panel);
            app.RecordNSpinner.Limits = [1 100];
            app.RecordNSpinner.ValueChangedFcn = createCallbackFcn(app, @RecordNSpinnerValueChanged, true);
            app.RecordNSpinner.FontName = 'Candara';
            app.RecordNSpinner.FontSize = 14;
            app.RecordNSpinner.Position = [98 364 70 22];
            app.RecordNSpinner.Value = 1;

            % Create DelayEditFieldLabel
            app.DelayEditFieldLabel = uilabel(app.Panel);
            app.DelayEditFieldLabel.HorizontalAlignment = 'right';
            app.DelayEditFieldLabel.FontName = 'Candara';
            app.DelayEditFieldLabel.FontSize = 14;
            app.DelayEditFieldLabel.Position = [29 187 55 22];
            app.DelayEditFieldLabel.Text = 'Delay(%)';

            % Create DelayEditField
            app.DelayEditField = uieditfield(app.Panel, 'numeric');
            app.DelayEditField.ValueChangedFcn = createCallbackFcn(app, @DelayEditFieldValueChanged, true);
            app.DelayEditField.FontName = 'Candara';
            app.DelayEditField.FontSize = 14;
            app.DelayEditField.Position = [98 187 67 22];
            app.DelayEditField.Value = 30;

            % Create ExpGainEditFieldLabel
            app.ExpGainEditFieldLabel = uilabel(app.Panel);
            app.ExpGainEditFieldLabel.HorizontalAlignment = 'right';
            app.ExpGainEditFieldLabel.FontName = 'Candara';
            app.ExpGainEditFieldLabel.FontSize = 14;
            app.ExpGainEditFieldLabel.Position = [30 158 54 22];
            app.ExpGainEditFieldLabel.Text = 'ExpGain';

            % Create ExpGainEditField
            app.ExpGainEditField = uieditfield(app.Panel, 'numeric');
            app.ExpGainEditField.ValueChangedFcn = createCallbackFcn(app, @ExpGainEditFieldValueChanged, true);
            app.ExpGainEditField.FontName = 'Candara';
            app.ExpGainEditField.FontSize = 14;
            app.ExpGainEditField.Position = [98 158 67 22];
            app.ExpGainEditField.Value = 5;

            % Create LCFHzEditFieldLabel
            app.LCFHzEditFieldLabel = uilabel(app.Panel);
            app.LCFHzEditFieldLabel.HorizontalAlignment = 'right';
            app.LCFHzEditFieldLabel.FontName = 'Candara';
            app.LCFHzEditFieldLabel.FontSize = 14;
            app.LCFHzEditFieldLabel.Position = [31 129 52 22];
            app.LCFHzEditFieldLabel.Text = 'LCF(Hz)';

            % Create LCFHzEditField
            app.LCFHzEditField = uieditfield(app.Panel, 'numeric');
            app.LCFHzEditField.ValueChangedFcn = createCallbackFcn(app, @LCFHzEditFieldValueChanged, true);
            app.LCFHzEditField.FontName = 'Candara';
            app.LCFHzEditField.FontSize = 14;
            app.LCFHzEditField.Position = [99 129 66 22];
            app.LCFHzEditField.Value = 250;

            % Create HCFHzEditFieldLabel
            app.HCFHzEditFieldLabel = uilabel(app.Panel);
            app.HCFHzEditFieldLabel.HorizontalAlignment = 'right';
            app.HCFHzEditFieldLabel.FontName = 'Candara';
            app.HCFHzEditFieldLabel.FontSize = 14;
            app.HCFHzEditFieldLabel.Position = [28 100 55 22];
            app.HCFHzEditFieldLabel.Text = 'HCF(Hz)';

            % Create HCFHzEditField
            app.HCFHzEditField = uieditfield(app.Panel, 'numeric');
            app.HCFHzEditField.ValueChangedFcn = createCallbackFcn(app, @HCFHzEditFieldValueChanged, true);
            app.HCFHzEditField.FontName = 'Candara';
            app.HCFHzEditField.FontSize = 14;
            app.HCFHzEditField.Position = [98 100 67 22];
            app.HCFHzEditField.Value = 2000;

            % Create FilterSwitchLabel
            app.FilterSwitchLabel = uilabel(app.Panel);
            app.FilterSwitchLabel.HorizontalAlignment = 'center';
            app.FilterSwitchLabel.FontName = 'Candara';
            app.FilterSwitchLabel.FontSize = 14;
            app.FilterSwitchLabel.Position = [36 245 41 22];
            app.FilterSwitchLabel.Text = 'Filter?';

            % Create FilterSwitch
            app.FilterSwitch = uiswitch(app.Panel, 'rocker');
            app.FilterSwitch.Items = {'No', 'Yes'};
            app.FilterSwitch.Orientation = 'horizontal';
            app.FilterSwitch.ValueChangedFcn = createCallbackFcn(app, @FilterSwitchValueChanged, true);
            app.FilterSwitch.FontName = 'Candara';
            app.FilterSwitch.FontSize = 14;
            app.FilterSwitch.Position = [106 249 29 13];
            app.FilterSwitch.Value = 'No';

            % Create UploadFilesButton
            app.UploadFilesButton = uibutton(app.Panel, 'push');
            app.UploadFilesButton.ButtonPushedFcn = createCallbackFcn(app, @UploadFilesButtonPushed, true);
            app.UploadFilesButton.IconAlignment = 'center';
            app.UploadFilesButton.BackgroundColor = [0.9412 0.9412 0.9412];
            app.UploadFilesButton.FontName = 'Candara';
            app.UploadFilesButton.FontSize = 14;
            app.UploadFilesButton.Position = [42 574 112 26];
            app.UploadFilesButton.Text = 'Upload Files';

            % Create PrintReportButton
            app.PrintReportButton = uibutton(app.Panel, 'push');
            app.PrintReportButton.ButtonPushedFcn = createCallbackFcn(app, @PrintReportButtonPushed, true);
            app.PrintReportButton.FontName = 'Candara';
            app.PrintReportButton.FontSize = 14;
            app.PrintReportButton.Position = [45 8 100 25];
            app.PrintReportButton.Text = 'Print Report';

            % Create DateDatePickerLabel
            app.DateDatePickerLabel = uilabel(app.Panel);
            app.DateDatePickerLabel.HorizontalAlignment = 'right';
            app.DateDatePickerLabel.FontName = 'Candara';
            app.DateDatePickerLabel.FontSize = 14;
            app.DateDatePickerLabel.Position = [18 40 34 22];
            app.DateDatePickerLabel.Text = 'Date';

            % Create DateDatePicker
            app.DateDatePicker = uidatepicker(app.Panel);
            app.DateDatePicker.FontName = 'Candara';
            app.DateDatePicker.FontSize = 14;
            app.DateDatePicker.Position = [67 42 115 22];
            app.DateDatePicker.Value = datetime([2023 6 9]);

            % Create Step1Label
            app.Step1Label = uilabel(app.Panel);
            app.Step1Label.BackgroundColor = [0.9412 0.9412 0.9412];
            app.Step1Label.FontName = 'Candara';
            app.Step1Label.FontSize = 14;
            app.Step1Label.Tooltip = {'Upload .txt files corresponding to the single pile test'};
            app.Step1Label.Position = [17 608 44 22];
            app.Step1Label.Text = 'Step  1';

            % Create Step2Label
            app.Step2Label = uilabel(app.Panel);
            app.Step2Label.BackgroundColor = [0.9412 0.9412 0.9412];
            app.Step2Label.FontName = 'Candara';
            app.Step2Label.FontSize = 14;
            app.Step2Label.Tooltip = {'Enter pile name, speed, frequency and length'};
            app.Step2Label.Position = [17 544 42 22];
            app.Step2Label.Text = 'Step 2';

            % Create Step3Label
            app.Step3Label = uilabel(app.Panel);
            app.Step3Label.BackgroundColor = [0.9412 0.9412 0.9412];
            app.Step3Label.FontName = 'Candara';
            app.Step3Label.FontSize = 14;
            app.Step3Label.Tooltip = {'Change test number'};
            app.Step3Label.Position = [18 394 42 22];
            app.Step3Label.Text = 'Step 3';

            % Create Step5Label_2
            app.Step5Label_2 = uilabel(app.Panel);
            app.Step5Label_2.BackgroundColor = [0.9412 0.9412 0.9412];
            app.Step5Label_2.FontName = 'Candara';
            app.Step5Label_2.FontSize = 14;
            app.Step5Label_2.Tooltip = {'Filter the stack'};
            app.Step5Label_2.Position = [17 274 43 22];
            app.Step5Label_2.Text = 'Step 5';

            % Create Step6Label_2
            app.Step6Label_2 = uilabel(app.Panel);
            app.Step6Label_2.BackgroundColor = [0.9412 0.9412 0.9412];
            app.Step6Label_2.FontName = 'Candara';
            app.Step6Label_2.FontSize = 14;
            app.Step6Label_2.Tooltip = {'Enter delay, gain, and frequency data for the filters'};
            app.Step6Label_2.Position = [15 216 43 22];
            app.Step6Label_2.Text = 'Step 6';

            % Create Step7Label
            app.Step7Label = uilabel(app.Panel);
            app.Step7Label.BackgroundColor = [0.9412 0.9412 0.9412];
            app.Step7Label.FontName = 'Candara';
            app.Step7Label.FontSize = 14;
            app.Step7Label.Tooltip = {'Modify date and export stack signal processing graph'};
            app.Step7Label.Position = [19 71 42 22];
            app.Step7Label.Text = 'Step 7';

            % Create Step4Label_2
            app.Step4Label_2 = uilabel(app.Panel);
            app.Step4Label_2.BackgroundColor = [0.9412 0.9412 0.9412];
            app.Step4Label_2.FontName = 'Candara';
            app.Step4Label_2.FontSize = 14;
            app.Step4Label_2.Tooltip = {'Choose records'};
            app.Step4Label_2.Position = [18 334 43 22];
            app.Step4Label_2.Text = 'Step 4';

            % Create StackButton
            app.StackButton = uibutton(app.Panel, 'push');
            app.StackButton.ButtonPushedFcn = createCallbackFcn(app, @StackButtonPushed, true);
            app.StackButton.FontName = 'Candara';
            app.StackButton.FontSize = 14;
            app.StackButton.Position = [48 301 100 25];
            app.StackButton.Text = 'Stack?';

            % Create MADEWITHMATLABSTUDENTLICENSEBYERIKASALAZARLabel
            app.MADEWITHMATLABSTUDENTLICENSEBYERIKASALAZARLabel = uilabel(app.UIFigure);
            app.MADEWITHMATLABSTUDENTLICENSEBYERIKASALAZARLabel.FontName = 'Candara';
            app.MADEWITHMATLABSTUDENTLICENSEBYERIKASALAZARLabel.FontSize = 9;
            app.MADEWITHMATLABSTUDENTLICENSEBYERIKASALAZARLabel.FontWeight = 'bold';
            app.MADEWITHMATLABSTUDENTLICENSEBYERIKASALAZARLabel.Position = [487 16 345 22];
            app.MADEWITHMATLABSTUDENTLICENSEBYERIKASALAZARLabel.Text = 'This application is made under a student license and its commercialization is prohibited';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = esappm

            runningApp = getRunningApp(app);

            % Check for running singleton app
            if isempty(runningApp)

                % Create UIFigure and components
                createComponents(app)

                % Register the app with App Designer
                registerApp(app, app.UIFigure)
            else

                % Focus the running singleton app
                figure(runningApp.UIFigure)

                app = runningApp;
            end

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end