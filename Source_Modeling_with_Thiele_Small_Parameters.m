% Abril 12 2020

% Este código calcula la Impedancia Electrica, Presion, Velocidad Volumetrica
% y de particula de N pistones empotrados en una pared infinita acoplados o no 
% a una linea de transmision a partir de sus Parametros Thiele Small a los cuales
% Es posible Inducir Variabilidad. 

clear;clc; close all;

% Número de pistones

N = 10;  
NN = 1:N;

% Thiele Small Dayton Audio PC68-4 (Se les Induce Variabilidad)

Blb = 1.93;                  % Parametro medido(media)
Bla = Blb*0.306;             % Blb*0.306; % Desviacion stnd
Bl = Bla.*randn(N,1) + Blb;  % [Tm] para no variar cambiar Bla por 0

Reb = 3.4;                   % Parametro medido(media)
Rea = 0;                     % Reb*0.0073; % Desviacion stnd
Re =  Rea.*randn(N,1) + Reb; % [Ohms] para no vari ar cambiar Rea por 0

Le  = (0.02e-3)*ones(N,1);   % [mH]

Mmsb = 0.0014;                  % Parametro medido(media)
Mmsa = 0;                       % Mmsb*0.0512; % Desviacion stnd
Mms  = Mmsa.*randn(N,1) + Mmsb; % [kg] para no variar cambiar Mmsa por 0

Cmcb = 0.9077e-3;               % Parametro medido(media)
Cmca = 0;                       %Cmcb*0.0757;% Desviacion stnd
Cmc = Cmca.*randn(N,1) + Cmcb;  % [mm/N] para no variar cambiar Cmca por 0

Rmcb = 0.4282;                  % Parametro medido(media)
Rmca = 0;                       % Rmcb*0.0747;% Desviacion stnd
Rmc = Rmca.*randn(N,1) + Rmcb;  % [Ns/m] para no variar cambiar Rmca por 0

Sd  = 0.00152*ones(N,1);        % [m^2]

Fs = 1./(2*pi*sqrt(Cmc.*Mms));  %[Hz]

% OCTAVAS

OCTAVAS = [31.5 63 125 250 500 1000 2000 4000 8000 16000];  % Octavas
fmin = round(OCTAVAS./sqrt(2));                             % Frecuencia inferior octavas
fmax = round(OCTAVAS.*sqrt(2));                             % Frecuencia superior octavas

% Ecuacion Que relaciona la Velocidad del Sonido con la Temperatura

Tc = 21.8255;                    % Temperatura [°C]
Tk = 273.15 + Tc;                % Temperatura Kelvin
Y  = 1.4;                        % Constante adiabática 
R  = 8.31;                       % Constante de Gas
M  = 29e-3;                      % Masa Molecular del Gas [kg/mol]
c  = sqrt((Y*R*Tk)/(M));         % Velocidad del Sonido
ro = 1.2041;                     % Densidad del Medio [kg/m^3]
Zcaracteristica = ro*c;          % Impedancia Caracteristica

% Creando vecores 

df = 1;                          % Periodo de Muestreo
f = min(fmin):df:max(fmax);      % Vector de frecuencias [HZ]  
w  = (2*pi).*f;                  % Vector de Frecuencia Angular [Rad/seg]
w = repmat(w,[N,1]);             % Matriz de frecuencia angular
k  = w./c;                       % Numero de Onda
l = 0.01;                        % Longitud del Tubo [m]

% Caracteristicas del Sistema

Ze = Re + (1i.*w.*Le);                                      % Impedancia electrica 
a = sqrt(Sd/pi);                                            % Radio del Altavoz [m]
ka = k.*a;                                                  % Radio del Piston x Numero de onda
aT = 0.045;                                                 % Radio del tubo [m] 
ST = pi*(aT^2);                                             % Area del tubo [m^2]
Mmd = Mms - ((16*ro*a.^3)/3);                               % Extraer la masa de aire que se acopla al diafragma
Zm =  Rmc + (1i.*w.*Mmd) + (1./(1i.*w.*Cmc));               % Impedancia Mecanica Total
Zam = (((((1/4)*((k*aT).^2))+(1i*0.61*k*aT))*(ro*c))/(ST)); % CONDICION DE FRONTERA
Zs = Zcaracteristica/ST;                                    % Impedancia caracteristica en una superficie determinada

% Matrices de Propagacion

n = (21*(1.84e-5))/Tc;                    % Coeficiente de viscosidad                    
alpha = (1/aT)*(sqrt((w*n)/(2*ro*c^2)));  % Coeficiente de atenuación
gama = k -(1i*alpha);                     % Coeficiente modelo disipativo 

% Coeficientes de Propagacion

A11 = ((Zam.*(Ze.*(Sd.^2).*Zs.*cos(gama*l) +(1i*sin(gama*l)).*(Zm.*Ze+(Bl.^2))) + Zs.*(1i*Zs.*Ze.*(Sd.^2).*sin(gama*l)+cos(gama*l).*(Ze.*Zm+(Bl.^2))))./(Zs.*Bl.*Zam.*Sd));
A12 = (1i*Zs.*Ze.*(Sd.^2).*sin(gama*l)+cos(gama*l).*((Bl.^2)+Ze.*Zm))./(Bl.*Sd);
A21 = ((Zam.*(Zs.*(Sd.^2).*cos(gama*l)+1i.*Zm.*sin(gama*l)) + Zs.*(1i*Zs.*(Sd.^2).*sin(gama*l)+Zm.*cos(gama*l)))./(Zs.*Bl.*Sd.*Zam));
A22 = (1i*Zs*(Sd.^2).*sin(gama*l)+Zm.*cos(gama*l))./(Bl.*Sd);

% ZE Y PZAM

Pref = 20e-6;                                       % Presion de Referencia
V = 2.83*ones(1,length(f));                         % Vector de voltaje constante
ZE = A11./A21;                                      % IMPEDANCIA ELECTRICA 
Phase_ZE = atan(imag(ZE)./real(ZE));                % FASE DE LA IMPEDANCIA ELECTRICA [Rad]
Phase_ZE_Deg = ((Phase_ZE*180)/pi);                 % FASE DE LA IMPEDANCIA ELECTRICA [Deg]
I = V./ZE;                                          % Corriente 
WE = I.*V;                                          % Potencia Electrica
Phase_WE = (atan(imag(WE)./real(WE))*180)/pi ;      % Fase de la potencia Electrica [Deg]

PZAM = V./A11;                                 % Precion en la Boca de la Linea de Transmicion
PZAMdB = 20*log10(abs(PZAM)/Pref);             % Precion en la Boca de la Linea de Transmicion en dB
Phase_PZAM = atan(imag(PZAM)./real(PZAM));     % Fase de la Precion en la Boca de la Linea de Transmicion [Rad]
Phase_PZAM_Deg = ((Phase_PZAM*180)/pi);        % Fase de la Precion en la Boca de la Linea de Transmicion [Deg]
U0 = (PZAM./(Sd.*Zam));                        % Velocidad de particula en el extremo abierto del tubo (Lossly Pipe) 
V_U0 = (PZAM./Zam);                            % Velocidad Volumentrica en el extremo abierto del tubo (Lossly Pipe) 
WA = PZAM.*V_U0;                               % Potencia Acustica
Phase_WA = (atan(imag(WA)./real(WA))*180)/pi;  % Fase de la Potencia Acustica 
Efficiency = WA./WE;                           % Relacion entre Potencia Electrica y Potencia Acustica 
Phase_Efficiency = (atan(imag(Efficiency)./...
    real(Efficiency))*180)/pi;                 % Fase de la Eficiencia 
X_U0V = V_U0./(1i*w(1,:));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;

Xplot=100;
Yplot=70;
width=1200;
height=600;

% FASE Y MAGNITUD DE LA IMPEDANCIA ELECTRICA DEL ALTAVOZ

figure(1)
subplot(211)
semilogx(f,abs(ZE),'Linewidth',2)
title('Phase and Magnitude Diagram Electrical Impedance','Fontname','Times','FontSize', 14)
ylabel('Electrical Impedance (\Omega)','Fontname','Times','FontSize', 13,'fontweight','bold')
ldg = legend(num2str(NN'));
title(ldg,'Source')
axis('tight')
grid on

subplot(212)
semilogx(f,Phase_ZE_Deg,'Linewidth',2)
xlabel('Frequency [Hz]','Fontname','Times','FontSize', 13,'fontweight','bold')
ylabel('Phase [Deg]','Fontname','Times','FontSize', 13,'fontweight','bold')
ldg = legend(num2str(NN'));
title(ldg,'Source')
axis('tight')
grid on

set(gcf,'position',[Xplot,Yplot,width,height])

% FASE Y MAGNITUD DE LA PRESION EN LA BOCA DE LA LINEA DE TRANSMICION

figure(2)
subplot(211)
semilogx(f,PZAMdB,'Linewidth',2)
title('Phase and Magnitude Diagram Mouth Pressure ','Fontname','Times','FontSize', 14)
ylabel(' SPL [dB] ','Fontname','Times','FontSize', 13,'fontweight','bold')
ldg = legend(num2str(NN'));
title(ldg,'Source')
axis('tight')
grid on

subplot(212)
semilogx(f,Phase_PZAM_Deg,'Linewidth',2)
xlabel('Frequency [Hz]','Fontname','Times','FontSize', 13,'fontweight','bold')
ylabel('Phase [Deg]','Fontname','Times','FontSize', 13,'fontweight','bold')
ldg = legend(num2str(NN'));
title(ldg,'Source')
axis('tight')
grid on

set(gcf,'position',[Xplot,Yplot,width,height])


% VARIABLES ELECTRICAS

figure(3)

subplot(311)
semilogx(f,abs(V),'Linewidth',2)
title('Electrical Variables','Fontname','Times','FontSize', 12,'fontweight','bold')
ylabel('Voltage [V]','Fontname','Times','FontSize', 13,'fontweight','bold')
legend('All Sources')
axis('tight'),grid on

subplot(312)
semilogx(f,abs(I),'Linewidth',2)
ldg = legend(num2str(NN'),'NumColumns',2);
title(ldg,'Source')
ylabel('Current [A]','Fontname','Times','FontSize', 13,'fontweight','bold')
axis('tight'),grid on

subplot(313)
semilogx(f,abs(ZE),'Linewidth',2)
ylabel('Impedance [\Omega]','Fontname','Times','FontSize', 13,'fontweight','bold')
ldg = legend(num2str(NN'),'NumColumns',2);
title(ldg,'Source')
axis('tight'),grid on

set(gcf,'position',[Xplot,Yplot,width,height])

% Eficiecia 

figure(4)

subplot(321)
semilogx(f,abs(WE),'Linewidth',2)
title('Magnitude','Fontname','Times','FontSize', 12,'fontweight','bold')
ylabel('Electric Power [W]','Fontname','Times','FontSize', 13,'fontweight','bold')
ldg = legend(num2str(NN'),'NumColumns',2);
title(ldg,'Source')
axis('tight'),grid on

subplot(322)
semilogx(f,Phase_WE,'Linewidth',2)
title('Phase [Deg]','Fontname','Times','FontSize', 12,'fontweight','bold')
ldg = legend(num2str(NN'),'NumColumns',2);
title(ldg,'Source')
axis('tight'),grid on

subplot(323)
semilogx(f,abs(WA),'Linewidth',2)
ylabel('Acoustic Power [W]','Fontname','Times','FontSize', 13,'fontweight','bold')
ldg = legend(num2str(NN'),'NumColumns',2);
title(ldg,'Source')
axis('tight'),grid on

subplot(324)
semilogx(f,Phase_WA,'Linewidth',2)
ldg = legend(num2str(NN'),'NumColumns',2);
title(ldg,'Source')
axis('tight'),grid on

subplot(325)
semilogx(f,abs(Efficiency),'Linewidth',2),axis tight, grid on
xlabel('Frequency [Hz]','Fontname','Times','FontSize', 13,'fontweight','bold')
ylabel('Efficiency','Fontname','Times','FontSize', 13,'fontweight','bold')
ldg = legend(num2str(NN'),'NumColumns',2);
title(ldg,'Source')
axis('tight'),grid on

subplot(326)
semilogx(f,Phase_Efficiency,'Linewidth',2),axis tight, grid on
xlabel('Frequency [Hz]','Fontname','Times','FontSize', 13,'fontweight','bold')
ldg = legend(num2str(NN'),'NumColumns',2);
title(ldg,'Source')
axis('tight'),grid on

set(gcf,'position',[Xplot,Yplot,width,height])


figure(5)

subplot(211)
semilogx(f,abs(V_U0),'Linewidth',2)
title('Spectrum Magnitude','Fontname','Times','FontSize', 12,'fontweight','bold')
ylabel('Volume velocity [m^3/s]','Fontname','Times','FontSize', 13,'fontweight','bold')
ldg = legend(num2str(NN'),'NumColumns',2);
title(ldg,'Source')
axis tight,grid on

subplot(212)
semilogx(f,abs(X_U0V),'Linewidth',2)
ylabel('Displacement [m]','Fontname','Times','FontSize', 13,'fontweight','bold')
xlabel('Frequency [Hz]','Fontname','Times','FontSize', 13,'fontweight','bold')
ldg = legend(num2str(NN'),'NumColumns',2);
title(ldg,'Source')
axis tight,grid on

set(gcf,'position',[Xplot,Yplot,width,height])