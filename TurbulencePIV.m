clear all, clc, clear memory
load('H:\Paper\Congresos\H1_3-7P1_Latam.mat')
[r, s]=size(u_original);

rho=1000;
mu=1.002*10^-3; % m2/s
X=x{1,1}(1,:);
yh=y{1,1}(:,1);
Yo=(max(yh)-yh);
ds=0.0010; %diametro de las gravas metros
Y=(Yo-ds)/max(Yo-ds);

for l=1:r
e(:,:,l)= u_original{l,1}-u_filtered{l,1}; %e matriz para determinar porcentaje de datos erroneos
end
    [a, b,c]=size(e);


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
       
p=vect>0.7;
for i=1:length(p)
    if p(i)==1
    u(:,:,i)=u_filtered{i,1};
    v(:,:,i)=v_filtered{i,1};

end
end
[a1 a2 a3]=size(u)
    for k=1:a3
               U(:,:,k)=u(:,:,k).*mask(:,:);
               V(:,:,k)=v(:,:,k).*mask(:,:);

        if U(:,:,k)==0
            U(:,:,k)=NaN;
        end
        if V(:,:,k)==0
            V(:,:,k)=NaN;
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
%     figure
%     plot(Uf(:,:,k),Vf(:,:,k),'*')
%     hold on
end


for k=1:a3
    UfVf(:,:,k)=Uf(:,:,k).*Vf(:,:,k);
    Uf2(:,:,k)=Uf(:,:,k).^2;
    Vf2(:,:,k)=Vf(:,:,k).^2;
end

%% Tildas
for i=1:a
    for j=1:b
Utilda(i,j)=Umed(i,j)-Upx(i,1);
Vtilda(i,j)=Vmed(i,j)-Vpx(i,1);
    end
end
UtVt=Utilda.*Vtilda;  
%% std tildas (esfuerzos de corte inducidos)

%% Intensidades Turbulentas y tau turbulento temporal

for i=1:a
    for j=1:b
rmsU(i,j)=(nanmean(Uf2(i,j,:))).^0.5; %Intensidad turbulenta U
rmsV(i,j)=(nanmean(Vf2(i,j,:))).^0.5; %Intensidad turbulenta V
ttt(i,j)=nanmean(UfVf(i,j,:)); %Tau turbulento promedio temporal

    end
end
for i=1:a
    rmsUp(i)=nanmean(rmsU(i,:));
    rmsVp(i)=nanmean(rmsV(i,:));
    tttp(i)=nanmean(ttt(i,:));
    Ftildas(i)=nanmean(UtVt(i,:));
    stdFtilda(i)=nanstd(UtVt(i,:));
end

%% Esfuerzos de corte
% Turbulento
for i=1:a
tt2p(1,i)= -rho*nanmean(ttt(i,:))%/(rho*9.81*R*ds/1000); %Tau turbulento doble promedio adm
ttildas(1,i)=-rho*Ftildas(i);
tt2p_std(1,i)=-rho*nanstd(ttt(i,:));
ttildas_std(1,i)=-rho*stdFtilda(i);
end

for i=1:length(tt2p)-1
    dtt(1)=0;
    dtt(i+1)=tt2p(i+1)-tt2p(i);
    dy(1)=0;
    dy(i+1)=Y(i+1)-Y(i);
end
m=dy./dtt;  
% if m(i)>100
%         m(i)=NaN
%     end
% Viscoso
for i=1:length(tt2p)-1
    dy(1,1)=0;
    du(1,1)=0;
    du(1,i+1)=-(Upx(i+1)-Upx(i));
    dy(1,i+1)=Y(i+1)-Y(i);
end

for i=1:a
    for j=1:b-1
    dys(1,j)=0;
    dus(1,j)=0;
    dus(i,j+1)=-(Umed(i,j+1)-Umed(i,j));
    end
    dustd(1,i)=nanstd(dus(i,:))
end

tvs=mu*dustd./dy;
tv=mu*du./dy;
ttotal=tt2p+ttildas+tv; %NO adimensional
ttotalstd=tt2p_std+ttildas_std+tvs;
%% Guardar datos

 save('H1-3-7_nb','p', 'X','Y','Uf','Vf','UfVf','Umed','Vmed','Upx','Vpx','ttt','tt2p','rmsU','rmsV', 'rmsUp','rmsVp','tttp','U', 'V', 'tv', 'ttotal','ttildas','ttotalstd','Yo')

%% Graficas
figure
subplot(2,2,1)
plot(ttildas,Y)
xlabel('u_{tilda}v_{tilda}')
ylabel('z/H')
grid on

subplot(2,2,2)
plot(Y,Upx)
xlabel('h/Y')
ylabel('Upx')
grid on

subplot(2,2,3)
plot(ttotal,Y)
hold on
% plot(ttotalstd+ttotal,Y)
% plot(-ttotalstd+ttotal,Y)
xlim([0,2])
xlabel('\tau_{total}')
ylabel('z/H')
grid on

% subplot(2,2,4)
% plot(tv,Y)
% ylabel('z/H')
% xlabel('\tau_v')
% grid on

% savefig('NOV18_P2acc.fig')
    
