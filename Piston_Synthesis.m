% Abril 12, 2020

% Este codigo sintetiza el campo radiado de N pistones con diferentes 
% parametros Thiele Small a partir de los calculos realizados en 
% Source_Modeling_with_Thiele_Small_Parameters.m

f_disp = 1000;
[~ , f_idx] = min(abs(f_disp - f));
H_dir = @(nu) 2*besselj(1,nu)./nu;  % función de directividad

res = 0.01;                         % Resolucion del dominio
dist_src = 0.089;                   % Separación entre pistones
kd = k.*dist_src;                   % Numero de onda por separacion entre fuentes
dim_dom = [6,6,3];                  % Dimensiones del dominio (m)
src_0 = [3,0,1.5];                  % Ubicación del piston (m)

% Distancia limite de los monopolos (Determinar que la longitud del arreglo no sobrepase las dimencione del dominio)
if ((N-1)*dist_src)/2 > dim_dom(1)
    disp (' El número de pistones excede las dimensiones del recinto o el monopolo se encuentra fuera del dominio');
    return
end

src = zeros(N,3); % Matriz de coordenadas de los monopolos
lentgh_array = (N)*dist_src;  % Longitud del arreglo
src(1:N,1) = (((-lentgh_array+dist_src)/2) : dist_src : (lentgh_array/2)); % Crea las cordenadas del arreglo para X ya que y es constante

% Sound Field reconstruction (CREANDO EL DOMINIO)

x_room = (0:res:dim_dom(1));
y_room = (0:res:dim_dom(2));
z_room = (0:res:dim_dom(3));
y_room = fliplr(y_room);

x_room = x_room -(src_0(1));
y_room = y_room-(src_0(2));
z_room = z_room-(src_0(3));
[X,Y] = meshgrid(x_room,y_room);

% Es necesario inicializar las diferentes magnitudes ya que cada uno de los pistones va a exibir una
% magnitud y un angulo diferente 

r = zeros(size(X,1),size(X,2),N);     % Inicializando matrices en 0
Synt_theta = zeros(size(X,1),size(X,2),N); % Inicializando matrices en 0
Np = zeros(length(x_room),length(y_room),N);    
NP = zeros(length(x_room),length(y_room));


tic
for n_idx = 1:size(src,1)
    
    x_scr = x_room - src(n_idx,1);
    y_scr = y_room - src(n_idx,2);

    [X_s,Y_s] = meshgrid(x_scr,y_scr);

    r = sqrt(X_s.^2 + Y_s.^2);
    Synt_theta = atan(X_s./Y_s);
    Synt_theta(isnan(Synt_theta)) = eps;
    Synt_theta(Synt_theta==0) = eps;

    Np(:,:,n_idx) =  U0(n_idx,f_idx).*(1j/2)*Zcaracteristica*a(n_idx,1).*ka(n_idx,f_idx).*...
                H_dir(ka(n_idx,f_idx).*sin(Synt_theta)).*...
                exp(-1j*k(n_idx,f_idx).*r)./r;   
                
    NP = NP + Np(:,:,n_idx);
                
end

N_Sourses_Radiation_dB = 20*log10(abs(NP)./Pref);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Xplot=5;
Yplot=20;
width=1400;
height=660;

figure(1)
subplot 121
pcolor(X,Y,real(NP));
shading interp;
title({'Synthesized Pressure Field (Real part)',['ka = ' num2str(ka(1,f_idx)')]...
    ,['kd = ' num2str(kd(1,f_idx)')]...
    ,['f = ',num2str(f(f_idx)),'Hz']}...
    ,'FontSize',17);
colorbar; set(gca, 'FontSize',21);
caxis([0 4]);
xlabel('X (m)','FontSize',24); ylabel('Y (m)','FontSize',24)
xlim([X(1,1,1), X(end,end,end)]); ylim([Y(end,end,end), Y(1,1,1)])
hold on;
scatter(src(:,1),src(:,2),'kd','filled')
colormap jet
subplot 122
pcolor(X,Y,N_Sourses_Radiation_dB);
title({'Synthesized SPL Field',['ka = ' num2str(ka(1,f_idx)')]...
    ,['kd = ',num2str(kd(1,f_idx)') ]...
    ,['f = ',num2str(f(f_idx)),'Hz']},'FontSize',17);
axis equal; caxis([80 115])
shading interp; colorbar; set(gca, 'FontSize',21);
xlabel('X (m)','FontSize',24); ylabel('Y (m)','FontSize',24)
xlim([X(1,1), X(end,end)]); ylim([Y(end,end), Y(1,1)])
hold on;
scatter(src(:,1),src(:,2),'kd','filled')
colormap jet
set(gcf,'position',[Xplot,Yplot,width,height])

figure(2)
subplot(121)
contourf(X,Y,real(NP));
shading interp; axis equal; 
title({'Synthesized Pressure Field -Real part-',['ka = ' num2str(ka(1,f_idx)')]...
    ,['kd = ',num2str(kd(1,f_idx)') ]...
    ,['f = ',num2str(f(f_idx)),'Hz']},'FontSize',17);
xlabel('X (m)','FontSize',18); ylabel('Y (m)','FontSize',18)
caxis([0 4]);
colormap jet ; colorbar;
hold on
scatter(src(:,1),src(:,2),'kd','filled')
hold off
set(gcf,'position',[Xplot,Yplot,width,height])

subplot(122)
contourf(X,Y,N_Sourses_Radiation_dB);
title({'Synthesized SPL Field',['ka = ' num2str(ka(1,f_idx)')]...
    ,['kd = ',num2str(kd(1,f_idx)') ]...
    ,['f = ',num2str(f(f_idx)),'Hz']},'FontSize',17);
xlabel('X (m)','FontSize',18); ylabel('Y (m)','FontSize',18)
shading interp; axis equal
caxis([80 115]); 
hold on
scatter(src(:,1),src(:,2),'kd','filled')
colormap jet; colorbar; 
hold off
set(gcf,'position',[Xplot,Yplot,width,height])