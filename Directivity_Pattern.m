% Abril 12, 2020

% Este código calcula el patron de directividad  por octavas de N pinstones
% empotrados en pantalla infinita a partir de la U0(Velocidad de particula) 
% estimada en Source_Modeling_with_Thiele_Small_Parameters.m, la sintesis 
% de sus patrones de directiviad, Su indice de directividad y su respuesta 
% en frecuencia total y por octavas a 0[°] y r[m]

% Sebastian Carvalho Salazar


H_dir = @(nu) 2*besselj(1,nu)./nu;  % función de directividad
r = 1;                              % distancia al observador (m)
d_theta = pi/180;                   % Diferencial de theta
theta = pi/2:-d_theta:-pi/2;        % Vector que representa el dominio polar (Semiesfera - Radiacion Frontal)
theta(isnan(theta)) = eps;          % Asignar a valores NaN el valor minimo capas de generar MATLAB
theta(theta == 0) = eps;            % Asignar a valores 0 el valor minimo capas de generar MATLAB
  
p = zeros(length(f),length(theta),N); % Patron de directividad de cada fuente
P = zeros(length(f),length(theta));   % Sintesis del campo radiado por las fuentes 
Q = zeros(1,length(f));               % indice de directividad
theta_disp = 0;                       % Theta en el que se quiere tomar la respuesta en frecuencia
[~ , theta_idx] = min(abs(theta_disp - theta)); % Triangular theta

OCT = zeros(length(OCTAVAS),N);                         % Presion en 0° por octavas[dB]     
P_OCT = zeros(length(OCTAVAS),length(theta),N);         % Patron de directividad por bandas de octava [dB]
KA_OCT = zeros(length(OCTAVAS),N);                      % promedio de todos los KA que componen cada octava
THETA = zeros(N,length(theta));
OCT_NSRC = zeros(length(OCTAVAS),1);
P_OCT_NSRC = zeros(length(OCTAVAS),length(theta));
Frequency_Response = zeros(N,length(f));

x = r.*sin(theta); % Cordenada Carteciana x
y = r.*cos(theta); % Cordenada Carteciana y
dist_src =    0.089;    % Distancia entre fuentes
lentgh_array = (N)*dist_src; % Longitud del arreglo
src(1:N,1) = (((-lentgh_array+dist_src)/2) : dist_src : (lentgh_array/2)); % Cordenadas del arreglo para X ya que Y es constante

tic
for N_SRC = 1:N
    
    for F = 1:df:length(f)
        
        x_src = x - src(N_SRC);
        y_src = y - 0;
       
        r = sqrt(x_src.^2+y_src.^2);
        theta = atan(x_src./y_src);
        theta(isnan(theta)) = eps;
        theta(theta==0) = eps;
        
        p(F,:,N_SRC) = U0(N_SRC,F).*(1j/2)*Zcaracteristica.*a(N_SRC,1).*ka(N_SRC,F).*...
                         H_dir(ka(N_SRC,F).*sin(theta)).*...
                            exp(-1j*k(N_SRC,F).*r)./r;  
                        
        Q(N_SRC,F) = (abs(p(F,theta_idx,N_SRC))^2) / (mean(abs(p(F,:,N_SRC)).^2));    
    
    end
    
    THETA(N_SRC,:) = theta;
    P = P + p(:,:,N_SRC);
    Frequency_Response(N_SRC,:) = 20.*log10(abs(p(:,theta_idx,N_SRC))./Pref);
    
    for i = 1:length(OCTAVAS)
        if   fmin(i) <= f(1,:) <= fmax(i)
            [~ , fmax_idx] = min(abs(fmax(i) - f));
            [~ , fmin_idx] = min(abs(fmin(i) - f));
            OCT(i,N_SRC) = mean(abs(p(fmin_idx:fmax_idx,theta_idx,N_SRC))); 
            P_OCT(i,:,N_SRC) = mean(abs(p(fmin_idx:fmax_idx,:,N_SRC)));
            KA_OCT(i,N_SRC) = mean(ka(N_SRC,fmin_idx:fmax_idx));
            
            OCT_NSRC(i) = mean(abs(P(fmin_idx:fmax_idx,theta_idx)));
            P_OCT_NSRC(i,:) = mean(abs(P(fmin_idx:fmax_idx,:)));
        end
        
    end
    
end
toc

DI = 10.*log10(Q);                                                   % Indice de directividad de cada fuente
OCT_NSRC_dB = 20*log10(OCT_NSRC./Pref);                              % Presion en octavas de la sintesis de N fuentes                     
OCT_dB = 20.*log10(OCT./Pref);                                       % Presion en 0° por octavas[dB]
P_OCT_dB = 20.*log10(P_OCT./Pref);                                   % Patron de directividad por bandas de octava [dB]
P_OCT_NSRC_dB = 20.*log10(P_OCT_NSRC./Pref);                         % Patron de directividad por bandas de octava de la sintesis de N fuentes
OCT_dB = [OCT_dB OCT_NSRC_dB];                                       % Une los SPL por octavas de cada fuente con el total 
Source_Radiation_dB = 20.*log10(abs(p)./Pref);                       % Presion de todas las fuentes [dB]
Frequency_Response_NSRC = 20.*log10(abs(P(:,theta_idx))./Pref);      % Presion en 0° para todas las frecuencias de la Sintesis de N fuentes
DI_NSRC = 10.*log10((abs(P(:,theta_idx)).^2)./mean(abs(P).^2,2));    % Indice de directividad de la sintesis de N fuentes
Source_Radiation_dB_T = 20.*log10(abs(P)./Pref);                     % Patron de directividad por frecuencia de cada fuente                 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Xplot=100;
Yplot=70;
width=1200;
height=600;

f_disp = 1000;
[~ , f_idx] = min(abs(f_disp - f));

figure(1) 
ax = polaraxes;
hold on
polarplot(mean(THETA,1),Source_Radiation_dB_T(f_idx,:),'k','Linewidth',2)
ax.ThetaZeroLocation = 'top';
thetalim(rad2deg([theta(end) theta(1)]))
title({'Directivity plot of the Sources Synthesis',...
    [' f = ',num2str(f(f_idx)),'Hz',...
    ' , ',...
    ' N = ', num2str(N),...
    ' , ',...
    'r  = 1','m',...
    ' , ',...
    'Distance Between Sources = ', num2str(dist_src),'m']})
thetaticks(-180:10:180)
rticks(0:10:max(max(P_OCT_NSRC_dB(:,:)))+10)
set(gcf,'position',[Xplot,Yplot,width,height])

figure(2)
semilogx(f,DI(:,:))
hold on
semilogx(f,DI_NSRC,'--','LineWidth',2)
title({'Directivity Index'...
        ,['Distance Between Sources = ', num2str(dist_src),'m']...
        ,['r = 1','m']...
        ,[' \theta = ',num2str(theta(theta_idx),2),'°']}...
        ,'FontName','Times New Roman','FontSize', 13)
xlabel('Frequency[Hz]','Fontname','Times','FontSize', 13,'fontweight','bold')
ylabel('Directivity [10log(Q)] ','Fontname','Times','FontSize', 13,'fontweight','bold')
ldg = legend([string(NN) 'Array'],'NumColumns',2,'Location','southwest');
title(ldg,'Source')
grid on,axis tight
hold off

set(gcf,'position',[Xplot,Yplot,width,height])

figure(3)
semilogx(f,Frequency_Response)
hold on
semilogx(f,Frequency_Response_NSRC,'--','LineWidth',2)
title({'Frequency Response'...
        ,['Distance Between Sources = ', num2str(dist_src),'m']...
        ,['r = 1','m']...
        ,[' \theta = ',num2str(theta(theta_idx),2),'°']}...
        ,'FontName','Times New Roman','FontSize', 13)
xlabel('Frequency[Hz]','Fontname','Times','FontSize', 13,'fontweight','bold')
ylabel('SPL[dB]','Fontname','Times','FontSize', 13,'fontweight','bold')
ldg = legend([string(NN) 'Array'],'NumColumns',2,'Location','southwest');
title(ldg,'Source')
grid on, axis tight
ylim([0 max(Frequency_Response_NSRC)+10])
hold off

set(gcf,'position',[Xplot,Yplot,width,height])

figure(4)
bar(OCT_dB)
title({'SPL Spectrum in Octaves '...
    ,['Distance Between Sources = ', num2str(dist_src),'m']...
    ,['r  = 1','m']...
    ,[' \theta = ',num2str(theta(theta_idx),2),'°']}...
    ,'FontName','Times New Roman','FontSize', 13)
xlabel('Octaves[Hz]','FontName','Times New Roman','FontSize', 12,'fontweight','bold')
ylabel('SPL[dB]','FontName','Times New Roman','FontSize', 11,'fontweight','bold')
xticks(1:length(OCTAVAS))
xticklabels(num2str(OCTAVAS'))
ldg = legend([string(NN) 'Array'],'NumColumns',2,'Location','northwest');
title(ldg,'Source')

set(gcf,'position',[Xplot,Yplot,width,height])

figure(5)

ax = polaraxes;
polarplot(mean(THETA),P_OCT_NSRC_dB(:,:))
ax.ThetaZeroLocation = 'top';
thetalim(rad2deg([theta(end) theta(1)]))
title({'Directivity pattern in octaves of the Sources Sythesis',...
[' N = ', num2str(N),...
' , ',...
'r  = 1','m',...
' , ',...
'Distance Between Sources = ', num2str(dist_src),'m']})
 ldg = legend(num2str(OCTAVAS',5),'NumColumns',1);
title(ldg,'Octaves')
set(gcf,'position',[Xplot,Yplot,width,height])
thetaticks(-180:10:180)
rticks(0:10:max(max(P_OCT_NSRC_dB(:,:)))+10)

figure(6)

Xplot=100;
Yplot=70;
width=1200;
height=600;

P_THETA = mean(THETA);

subplot(121)
contourf((P_THETA*180)./pi,f,abs(P))
title({'Pressure Magnitude in PA'...
    ,['Distance Between Sources = ', num2str(dist_src),'m']...
    ,['r = 1','m']}...
    ,'FontName','Times New Roman','FontSize', 13)
shading interp;
annotation('textbox',[.9 .5 .1 .2],'String','Text outside the axes','EdgeColor','none')
ylabel('Frequency[Hz]','Fontname','Times','FontSize', 13,'fontweight','bold')
xlabel('Angle [Deg]','Fontname','Times','FontSize', 13,'fontweight','bold')
caxis([2e-6 5])
colormap jet; colorbar

subplot(122)
contourf((mean(THETA)*180)./pi,f,Source_Radiation_dB_T)
title({'SPL Magnitude '...
    ,['Distance Between Sources = ', num2str(dist_src),'m']...
    ,['r = 1','m']}...
    ,'FontName','Times New Roman','FontSize', 13)
ylabel('Frequency[Hz]','Fontname','Times','FontSize', 13,'fontweight','bold')
xlabel('Angle [Deg]','Fontname','Times','FontSize', 13,'fontweight','bold')
shading interp;
caxis([20 100])
colormap jet; colorbar

set(gcf,'position',[Xplot,Yplot,width,height])