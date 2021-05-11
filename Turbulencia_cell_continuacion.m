%% Análisis PIV en cuartiles de las fluctuaciones turbulentas en tres puntos de la eyección
% Analizar por cuartiles los datos tomados con PIV, caracteriza
% elipses considerando: Como ángulo de inclinación la pendiente de la linea de tendencia de la serie de datos, como radio mayor y radio menor. 
%  
%% Procesamiento de resultados PIV
clear all, clc, clear memory
%H:\Natalia\Montaje2\matlab_NOV19pared_eyeccion.mat
load('H:\Paper\Congresos\piv_may5_centro_dcc.mat');
[r, s]=size(u_original);
fs=1000;%Hz
rho=1000;
mu=1.002*10^-3; % m2/s
X=x{1,1}(1,:);
yh=y{1,1}(:,1);
Yo=(max(yh)-yh);

ds=0/1000; %diametro de las gravas metros
Y=(Yo-ds)/max(Yo-ds);

for l=1:r
e(:,:,l)= u_original{l,1}-u_filtered{l,1}; %e matriz para determinar porcentaje de datos erroneos
end
    [a, b,c]=size(e);
% t=[1:200:2*17460]/1000;

    for k=1:r
        m(k)=0; n(k)=0;
        for i=1:a
            for j=1:b
                if e(i,j,k)==0
                    m(k)=m(k)+1;
                end
                    n(k)=sum(sum(isnan(e(:,:,k))));
            end
        end
        vect(k)=(m(k)+n(k))/(a*b);
    end

        for i=1:a
            for j=1:b
                if isnan(e(i,j,1));
                    mask(i,j)=0;
                else 
                    mask(i,j)=1;
                end
            end
        end
       
p=vect>0.8;
for i=1:length(p)
    if p(i)==1
    u(:,:,i)=u_filtered{i,1};
    v(:,:,i)=v_filtered{i,1};
    umag(:,:,i)=velocity_magnitude{i,1};
    end
end
[a1 a2 a3]=size(u)
    for k=1:a3
        U(:,:,k)=-u(:,:,k).*mask(:,:);
        V(:,:,k)=-v(:,:,k).*mask(:,:);
        Ucomp(:,:,k)=umag(:,:,k).*mask(:,:);
               
        if U(:,:,k)==0
            U(:,:,k)=NaN;
        end
        if V(:,:,k)==0
            V(:,:,k)=NaN;
        end
        if Ucomp(:,:,k)==0
            Ucomp(:,:,k)=NaN;
        end
    end
    
%% Velocidades medias
    for i=1:a
        for j=1:b
            Umed(i,j)=nanmean(U(i,j,:));
            Vmed(i,j)=nanmean(V(i,j,:));
        end
    end
%% Doble promedios
for i=1:a
    Upx(i,1)=nanmean(Umed(i,:));%Udoblepromedio
    Vpx(i,1)=nanmean(Vmed(i,:));%Vdoblepromedio
end    
%% Fluctuaciones
for k=1:a3
    Uf(:,:,k)=U(:,:,k)-Umed(:,:);
    Vf(:,:,k)=V(:,:,k)-Vmed(:,:);

end


for k=1:a3
    UfVf(:,:,k)=Uf(:,:,k).*Vf(:,:,k);
    Uf2(:,:,k)=Uf(:,:,k).^2;
    Vf2(:,:,k)=Vf(:,:,k).^2;
end


% figure
% for i=1:100:c
%     plot(Uf(:,:,i)/mean(mean(Umed)),Vf(:,:,i)/mean(mean(Umed)),'.b')
%     xlim([-1.5 1.5])
%     ylim([-1.5 1.5])
%     hold on
%     grid on
%     p = calculateEllipse(0, 0, radM/mean(mean(Umed)), radm/mean(mean(Umed)), ang, 100); %   p = calculateEllipse(0, 0, 0.023, 0.016, 7, 100); para medio poroso
%     plot(p(:,1), p(:,2), '.-g'), axis equal
%     xlabel('u^{''}')
%     ylabel('w^{''}')
%     title('NOV19PARED')
%     set(gca, 'fontsize', 14)
% end


%%  Hipotesis de turbulencia congelada de Taylor
imagesc(Uf(:,:,2728))

P=[1 2 3]; % Puntos
% coordenadas
%px=[3 10 17]; py =[18] Jorge
%px=[20 17 14]; py=[36] Nov19_lecho
%px=[21 16 11]; py=[15] Nov19_eyeccion
%px=[15 10 5]; py=[15]  Oct07
%px=[20 15 10]; py=[40]  Nov19 full water
for y=1:length(P)
    y=P(y)
   if y==1 %Aguas Arriba
       px=17; %jorge 3
       py=30; %40 water %jorge 18
   else if y==2 %Centro
       px=12; %jorge 10
       py=30; %40 water %jorge 18
       else if y==3  %Aguas abajo
       px=7; %jorge 17
       py=30; %10 water %jorge 18       
       end
       end
   end
m1=1;n1=1; m2=1;n2=1;
  for k=1:a3
      % Hipotesis de turbulencia congelada de Taylor
       Umagt(:,k)=Ucomp(:,px,k);
       Ut(:,k)=Uf(:,px,k); 
       Vt(:,k)=Vf(:,px,k);
       Vel_l(:,k)=U(:,px,k); %Magnitud de velocidad en la longitudinal
       Vel_y(:,k)=V(:,px,k); %Magnitud de velocidad en la vertical
       tadm(1,k)=k*nanmean(U(:,px,k))/(max(Y)); %tiempo adimensional
      
       U1_p(1,k)=U(py,px,k); %Serie de velocidad en un punto
       V1_p(1,k)=V(py,px,k); %Serie de velocidad en un punto
       U_p=U1_p(~isnan(U1_p));
       V_p=V1_p(~isnan(V1_p));       
       
       Q1f_u(1,k)=Uf(py,px,k);
       Q1f_v(1,k)=Vf(py,px,k);
       Qf_u=Q1f_u(~isnan(Q1f_u));
       Qf_v=Q1f_v(~isnan(Q1f_v));
       QUfQVf= Q1f_u.*Q1f_v;
%        [s b]=polyfit(Qf_u, Qf_v,1)
%        ang=-atan(s(1))*180/pi
       
 if UfVf(py,px,k)<0 & Uf(py,px,k)<0
        Q2(1,m1)=UfVf(py,px,k);
        Q2_K(1,m1)=UfVf(py,px,k)/nanmean(UfVf(py,px,:));
        Q2_u(1,m1)=Uf(py,px,k);
        Q2_v(1,m1)=Vf(py,px,k);
        m1=m1+1;
end
if UfVf(py,px,k)<0 & Uf(py,px,k)>0
        Q4(1,n1)=UfVf(py,px,k);
        Q4_K(1,n1)=UfVf(py,px,k)/nanmean(UfVf(py,px,:));
        Q4_u(1,n1)=Uf(py,px,k);
        Q4_v(1,n1)=Vf(py,px,k);
        n1=n1+1;
end
if UfVf(py,px,k)>0 & Uf(py,px,k)>0
        Q1(1,m2)=UfVf(py,px,k);
        Q1_K(1,m2)=UfVf(py,px,k)/nanmean(UfVf(py,px,:));
        Q1_u(1,m2)=Uf(py,px,k);
        Q1_v(1,m2)=Vf(py,px,k);
        m2=m2+1;
end
if UfVf(py,px,k)>0 & Uf(py,px,k)<0
        Q3(1,n2)=UfVf(py,px,k);
        Q3_K(1,n2)=UfVf(py,px,k)/nanmean(UfVf(py,px,:));
        Q3_u(1,n2)=Uf(py,px,k);
        Q3_v(1,n2)=Vf(py,px,k);
        n2=n2+1;
end   
  end
       [s b]=polyfit(Qf_u, Qf_v,1)
       ang=-atan(s(1))*180/pi
       
% Ordena en orden descendente las fluctuaciones turbulentas por cuartil.  
Q1s=sort(Q1_K,'descend');
Q2s=sort(Q2_K,'descend');
Q3s=sort(Q3_K,'descend');
Q4s=sort(Q4_K,'descend');

UfVfm=nanmean(nanmean(nanmean(UfVf)));
UfVfK=-QUfQVf/abs(UfVfm);
UfVfs=sort(UfVfK,'descend');


% Probabilidad

for i=1:length(Q1)
Q1p(i)=i/(1+length(Q1))*100;
end
for i=1:length(Q2)
Q2p(i)=i/(1+length(Q2))*100;
end
for i=1:length(Q3)
Q3p(i)=i/(1+length(Q3))*100;
end
for i=1:length(Q4)
Q4p(i)=i/(1+length(Q4))*100;
end
for i=1:length(UfVf)
UfVfp(i)=i/(1+length(UfVf))*100;
end

% Calcula elipse
F_Q1m=floor(Q1p)';
F_Q2m=floor(Q2p)';
F_Q3m=floor(Q3p)';
F_Q4m=floor(Q4p)';
F_QUfVfpm=floor(UfVfp)';

Q1m=Q1s'; Q1m(:,2)=Q1p'; Q1m(:,3)=F_Q1m; %Q1m(:,4)=Vfx1; Q1m(:,5)=Vfz1; 
Q2m=Q2s'; Q2m(:,2)=Q2p'; Q2m(:,3)=F_Q2m; %Q2m(:,4)=Vfx2; Q2m(:,5)=Vfz2;
Q3m=Q3s'; Q3m(:,2)=Q3p'; Q3m(:,3)=F_Q3m; %Q3m(:,4)=Vfx3; Q3m(:,5)=Vfz3;
Q4m=Q4s'; Q4m(:,2)=Q4p'; Q4m(:,3)=F_Q4m; %Q4m(:,4)=Vfx4; Q4m(:,5)=Vfz4;
QUfVfm(:,1)=UfVfs'; QUfVfm(:,2)=UfVfp; QUfVfm(:,3)=F_QUfVfpm;

Q1_95(:,1)=min(find(Q1m(:,3)==90));
Q2_5(:,1)=min(find(Q2m(:,3)==10));
Q3_95(:,1)=min(find(Q3m(:,3)==90));
Q4_5(:,1)=min(find(Q4m(:,3)==10));
QUfVfm_95(:,1)=min(find(QUfVfm(:,3)==90));
QUfVfm_5(:,1)=min(find(QUfVfm(:,3)==10));

% Dimension para Rb
Q1_K95=Q1m(Q1_95,1);
d1= min(find(Q1_K==Q1_K95));
Ru1b=Q1_u(d1);
Rv1b=Q1_v(d1);
Rb1=(Ru1b^2+Rv1b^2)^0.5

Q_K95m=QUfVfm(QUfVfm_95,1);
d1m=min(find(UfVfK==Q_K95m))
Ru1bm=Q1f_u(d1m);
Rv1bm=Q1f_v(d1m);
Rb1m=(Ru1bm^2+Rv1bm^2)^0.5

Q3_K95=Q3m(Q3_95,1);
d3= find(Q3_K==Q3_K95);
Ru3b=Q3_u(d3);
Rv3b=Q3_v(d3);
Rb3=(Ru3b^2+Rv3b^2)^0.5
Rb=min(Rb1,Rb3);

% Dimensiones para los diametros de las elipses - Ra
Q2_K5=Q2m(Q2_5,1);
d2= max(find(Q2_K==Q2_K5));
Ru2a=Q2_u(d2);
Rv2a=Q2_v(d2);
Ra2=(Ru2a^2+Rv2a^2)^0.5

Q4_K5=Q4m(Q4_5,1);
d4= max(find(Q4_K==Q4_K5));
Ru4a=Q4_u(d4);
Rv4a=Q4_v(d4);
Ra4=(Ru4a^2+Rv4a^2)^0.5
Ra=min(Ra2,Ra4);

Q_K5m=QUfVfm(QUfVfm_5,1);
d2m=min(find(UfVfK==Q_K5m))
Ru2bm=Q1f_u(d2m);
Rv2bm=Q1f_v(d2m);
Rb2m=(Ru2bm^2+Rv2bm^2)^0.5

%Radios
radM=max(Rb1m,Rb2m);
radm=min(Rb1m,Rb2m);
p = calculateEllipse(0, 0, radM , radm, ang, 100)


name=strcat('May05dcc_P',num2str(y))
save(name,'Yo','X','Vel_l','Vel_y','Qf_u','Qf_v','Q1','Q1_K','Q1_u','Q1_v','Q2','Q2_K','Q2_u','Q2_v','Q3','Q3_K','Q3_u','Q3_v','Q4','Q4_K','Q4_u','Q4_v','radM','radm','ang', 'Umagt','Ut','Vt','U','V','Ucomp','U_p','V_p')

% fileID = fopen(name,'w');
% fprintf(fileID,'%0s %10s %12s\n','Ra','Rb','Ángulo');
% fprintf(fileID,'%0.4f %9.4f %9.4f\n',radM, radm, ang);
% fclose(fileID);


FFT=fft2(Qf_u);
absFFT=abs(FFT).^2;
f=[1:1:length(Qf_u)]/fs;
f_1=f(25:length(Qf_u));

p1=polyfit(log(f), log(absFFT),1)
z=@(f_1) (10^p1(2))*f_1.^p1(1);
% fplot(z,[f_1(1),f_1(end)])
% fprintf('exponente a= %2.3f\n',p1(1));
% fprintf('coeficiente c = %3.3f\n',(10^p1(2)));

% figure
% loglog(f,absFFT)
% hold on
% grid on
% fplot(z,[f(25),f(end)])
% xlabel('frecuency [Hz]')
% ylabel('|FFT(u{''})|^2')
% title(name)

figure
plot(QUfVfm(:,1),QUfVfm(:,2),'-')
grid on
hold on
plot([-40 80],[10 10],':k')
plot([-40 80],[90 90],':k')
% plot(Q1m(:,1),Q1m(:,2),'.b')
% plot(Q2m(:,1),Q2m(:,2),'.r')
% plot(Q3m(:,1),Q3m(:,2),'.g')
% plot(Q4m(:,1),Q4m(:,2),'.k')
set(gca, 'yDir','normal','fontsize', 18)
xlabel('K')
ylabel('%')

fig=figure
plot(Q2_u/nanmean(nanmean(Vel_l)),Q2_v/nanmean(nanmean(Vel_l)),'.','MarkerSize',2)
hold on 
grid on
plot(Q1_u/nanmean(nanmean(Vel_l)),Q1_v/nanmean(nanmean(Vel_l)),'.','MarkerSize',2)
plot(Q3_u/nanmean(nanmean(Vel_l)),Q3_v/nanmean(nanmean(Vel_l)),'.','MarkerSize',2)
plot(Q4_u/nanmean(nanmean(Vel_l)),Q4_v/nanmean(nanmean(Vel_l)),'.','MarkerSize',2) %,'MarkerEdgeColor',[1 0.8 0.2])
set(gca, 'yDir','normal','fontsize', 18)
% plot(Ru1b,Rv1b,'pk')
% plot(Ru3b,Rv3b,'pk')
% plot(Ru2a,Rv2a,'ok')
% plot(Ru4a,Rv4a,'ok')
plot(p(:,1)/nanmean(nanmean(Vel_l)), p(:,2)/nanmean(nanmean(Vel_l)), '.-k')%, axis equal
ylim([-1.5 1.5])
xlim([-1.5 1.5])
% xlim([-0.04 0.04])
% ylim([-0.04 0.04])
xlabel('u^{''}/<U>')
ylabel('w^{''}/<U>')
title(name)

% rut1=strcat('H:\Natalia\Montaje2\nov19\',name,'.png')
% saveas(fig,rut1)

end

break
%% Cuartiles

%% Figuras
quiver(Ut, Vt,'k')

figure
imagesc(U(:,:,2728))
figure
plot(Q2_p1u,Q2_p1v,'.r')
hold on 
grid on
plot(Q1_p1u,Q1_p1v,'.b')
plot(Q3_p1u,Q3_p1v,'.g')
plot(Q4_p1u,Q4_p1v,'.k')

figure
plot(Qf_p1u,Qf_p1v,'.r')
% Cuartiles
figure
subplot(1,3,1)
plot(Q2_p1u,Q2_p1v,'.r')
hold on 
grid on
plot(Q1_p1u,Q1_p1v,'.b')
plot(Q3_p1u,Q3_p1v,'.g')
plot(Q4_p1u,Q4_p1v,'.k')
xlim([-0.04 0.04])
ylim([-0.04 0.04])
xlabel('u^{''}')
ylabel('w^{''}')
title('Punto1')

subplot(1,3,2)
plot(Q2_p2u,Q2_p2v,'.r')
hold on 
grid on
plot(Q1_p2u,Q1_p2v,'.b')
plot(Q3_p2u,Q3_p2v,'.g')
plot(Q4_p2u,Q4_p2v,'.k')
xlim([-0.04 0.04])
ylim([-0.04 0.04])
xlabel('u^{''}')
ylabel('w^{''}')
title('Punto2')

subplot(1,3,3)
plot(Q2_p3u,Q2_p3v,'.r')
hold on 
grid on
plot(Q1_p3u,Q1_p3v,'.b')
plot(Q3_p3u,Q3_p3v,'.g')
plot(Q4_p3u,Q4_p3v,'.k')
xlim([-0.04 0.04])
ylim([-0.04 0.04])
xlabel('u^{''}')
ylabel('w^{''}')
title('Punto3')

