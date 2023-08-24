% step 4
% transitorio iniziale, materiale 2 strati + parte pulsata+ evaporazione
close all
clearvars
clc

% set grafici
set(0,'defaultlinelinewidth',1)
set(0,'defaultaxesfontsize',16)

%proprietà e dimensioni
qminmin=1e6;
qmin=5e6; %W/m2
qmax1=15e6; % carico termico suoerficiale massimo
qmax2=35e6; % carico termico suoerficiale massimo
qmax3=60e6; % carico termico suoerficiale massimo

deltac=2e-3; % dt tra la fine di una pulsazione e l'inizio della successiva
deltap=1e-3; % delta t di pulsazione

deltacps=1e-3; %m
deltaBNC=2e-2;  %m
T_water=300;  %°C
hh=1;     % W/m2/K

%Litio_SN
cv_Li=1.6/1e-6;  %J/k/m3
k_Li=@(T) 27+T./40; % W/m/k

%Molibdeno
cv_Mo=2.6/1e-6; %J/k/m3
k_Mo=@(T) 145-T./33; % W/m/k

%Bulk
cv_BNC=2.57/1e-6;  %J/k/m3
k_BNC=16;% W/m/k

%CPS
f=0.5;
cv_CPS=cv_Li*f+cv_Mo*(1-f);
k_CPS=@(T) k_Li(T).*f+k_Mo(T).*(1-f);

%%% evaporazione
Pv=@(T) 10.^(10.061-8023./(T+273.15));
Nav=6.02214076e23; 
MM_Li=6.941; % g/mol
mLi=6.941/Nav*1e-3; % kg
kb=1.38065e-23; % J/k

g2Li=@(T) 1.66*1.9*Pv(T)./sqrt(2*pi*mLi*kb.*(273.15+T));     %%%fattore correttivo = 1
m2Li=@(T) g2Li(T)./Nav;
q2Li=@(T) m2Li(T)*136e3;

%%%%% evaporazione Sn
Nav=6.02214076e23; 
MM_tin=118.71; % g/mol
msn=118.71*1.66054e-27; % kg
kb=1.38065e-23; % J/k
pv=@(T) 10.^(5.006+5.262+(-15332./(T+273.15)));

g2Sn=@(T) 1.66*0.9*pv(T)./sqrt(2*pi*msn*kb.*(273.15+T));
m2Sn=@(T) g2Sn(T)./Nav;
q2Sn=@(T) m2Sn(T)*296e3;


q2=@(T) 0.27*q2Li(T)+.73*q2Sn(T);
% discretizzazione
dx1=1e-5;
x1=(0:dx1:deltacps)';
n1=length(x1);

dx2=5e-5;
x2=(deltacps+dx2:dx2:deltacps+deltaBNC)';
n2=length(x2);

xx=[x1;x2];
nn=length(xx);

dt=1e-3;
time(1)=0;

jj=1;
T0=350*ones(nn,1); % Tguess
toll=1e-8;
tollt=1e-6;
errt=1;
Tp=T0;  % Termine noto equazione
Told=T0; % Temperatura all'istante precedente
Ts(jj)=Tp(1); % Temperatura sulla superficie del CPS
Tinterf(jj)=Tp(n1); % Temperatura interfaccia CPS-BNC
Tfluid(jj)=Tp(end); % Temperatura all'interfaccia acqua Bulk 
for  ii=1:1.025/dt
    idx=0;
    jj=jj+1;
    time(jj)=time(jj-1)+dt;
    errc=1;
    
    while errc>toll
        idx=idx+1;
        %assemblaggio matrice
        Tcps=T0(1:n1);
        TBNC=T0(n1+1:end);
        
        %%% cps
        Tcpsp=[Tcps(2:end);Tcps(end)];
        Tcpsm=[Tcps(1);Tcps(1:end-1)];
        kicps=k_CPS(Tcps);
        kipcps=k_CPS(Tcpsp);
        kimcps=k_CPS(Tcpsm);

    
        sub_cps=-dt*1/dx1^2/cv_CPS*(kicps-(kipcps-kimcps)./4);
        main_cps=1+2/dx1^2*dt*kicps./cv_CPS;
        sup_cps=-dt*1/dx1^2/cv_CPS*(kicps+(kipcps-kimcps)./4);
        
        %%% Bulk
        aa=dt*k_BNC/dx2^2/cv_BNC;
        sub_b=-aa*ones(length(TBNC),1);
        main_b=1+2*aa*ones(length(TBNC),1);
        sup_b=sub_b;
        
        %%% dominio intero
        sub=[sub_cps(2:end);sub_b;0];
        main=[main_cps;main_b];
        sup=[0;sup_cps;sup_b(1:end-1)];
    
    
        AA=spdiags([sub,main,sup],-1:1,nn,nn);
        

        qev(jj)=q2(T0(1));
        %BC
        AA([1 n1 end],:)=0;
        
        AA(1,1:3)=[-3 4 -1];
        Tp(1)=-(qminmin-qev(jj))*2*dx1/k_CPS(T0(1));
        
        AA(n1,n1-1)=-kicps(end)/k_BNC*dx2/dx1;
        AA(n1,n1)=kicps(end)/k_BNC*dx2/dx1+1;
        AA(n1,n1+1)=-1;
        Tp(n1)=0;
    
        AA(end,end)=1;
        Tp(end)=330;
    
        TT=AA\Tp;
        errc=norm(TT-T0)/norm(T0-T_water);
        T0=TT;
    end
    Ts(jj)=TT(1);
    Tinterf(jj)=TT(n1);
    Tfluid(jj)=TT(end);
    errt=norm(TT-Told)/norm(Told-T_water);
    Tp=TT;
    Told=TT;
end




t0=time(end);
fun=@(time) abs(sin(pi/2*((time-t0)/5e-3)));
T0=TT;
T_b4_elms=TT;

figure(1)
plot(xx*1e3,TT,'displayname','Profilo prima del carico pulsato',LineWidth=1.8)
grid on
hold on
xlabel('spessore [mm]')
ylabel('Temperatura [°C]')


%%
dt=.5e-3;
idx=1;
ooo=1;
ppp=1;
istanti=[0.05 0.1 0.151];
for i=1:0.2/dt
    time(jj+i)=time(jj+i-1)+dt;
    if dt*ooo<=deltac
        qs(i)=qmin;
        ooo=ooo+1;
    else
        if dt*ooo-deltac<deltap
            if dt*i<=0.05 || time(jj+i)>=1.2
                qs(i)=qmax1;
            elseif dt*i<=0.1
                qs(i)=qmax2;
            elseif dt*i<0.2
                qs(i)=qmax3;
            end
            ooo=ooo+1;
        else
            qs(i)=qmin;
            ooo=1;
        end
    end
    errc=1;
    
    while errc>toll
        idx=idx+1;
       %assemblaggio matrice
        Tcps=T0(1:n1);
        TBNC=T0(n1+1:end);
        
        %%% cps
        Tcpsp=[Tcps(2:end);Tcps(end)];
        Tcpsm=[Tcps(1);Tcps(1:end-1)];
        kicps=k_CPS(Tcps);
        kipcps=k_CPS(Tcpsp);
        kimcps=k_CPS(Tcpsm);

    
        sub_cps=-dt*1/dx1^2/cv_CPS*(kicps-(kipcps-kimcps)./4);
        main_cps=1+2/dx1^2*dt*kicps./cv_CPS;
        sup_cps=-dt*1/dx1^2/cv_CPS*(kicps+(kipcps-kimcps)./4);
        
        %%% Bulk
        aa=dt*k_BNC/dx2^2/cv_BNC;
        sub_b=-aa*ones(length(TBNC),1);
        main_b=1+2*aa*ones(length(TBNC),1);
        sup_b=sub_b;
        
        %%% dominio intero
        sub=[sub_cps(2:end);sub_b;0];
        main=[main_cps;main_b];
        sup=[0;sup_cps;sup_b(1:end-1)];
    
    
        AA=spdiags([sub,main,sup],-1:1,nn,nn);
        

        qev(jj+i)=q2(T0(1));
        %BC
        AA([1 n1 end],:)=0;
        
        AA(1,1:3)=[-3 4 -1];
        Tp(1)=-(qs(i)-qev(jj+i))*2*dx1/k_CPS(T0(1));
        
        AA(n1,n1-1)=-kicps(end)/k_BNC*dx2/dx1;
        AA(n1,n1)=kicps(end)/k_BNC*dx2/dx1+1;
        AA(n1,n1+1)=-1;
        Tp(n1)=0;

        AA(end,end)=1;
        Tp(end)=330;
    
        TT=AA\Tp;
        errc=norm(TT-T0)/norm(T0-T_water);
        T0=TT;
    end
    if i*dt==istanti(ppp)
        figure(1)
       
        plot(xx*1e3,TT,'DisplayName',strcat('Profilo dopo ',num2str(i*dt),' s dall''inizio carico pulsato (VS)'),LineWidth=1.8)
        hold on
        if ppp<3
            ppp=ppp+1;
        else
            ppp=3;
        end
    end
    Ts(jj+i)=TT(1);
    Tinterf(jj+i)=TT(n1);
    Tfluid(jj+i)=TT(end);
    errt=norm(TT-Told)/norm(Told-T_water);
    Tp=TT;
    Told=TT;
end



figure(1)
plot(xx*1e3,TT,'DisplayName','Profilo dopo carico pulsato (VS)',LineWidth=1.8)
legend()

figure(2)
plot(time,Ts,'LineWidth',1,'DisplayName','Temperatura superficie CPS')
hold on
%plot(time,Tfluid,'LineWidth',1,'DisplayName','Temperatura interfaccia con accqua')
plot(time,Tinterf,'LineWidth',1,'DisplayName','Temperatura interfaccia CPS-W')
grid on
xlabel('tempo [s]')
ylabel('Temperatura [°C]')
legend('Location','best')


figure(3)
plot(time(jj+1:end)-t0,qs./1e6,'LineWidth',0.5)
grid on
title('Carico dovuto a ELMs')
xlabel('tempo [s]')
ylabel('q'''' [MW/m^2]')

figure(4)
plot(time(jj-10:end),Ts(jj-10:end),'LineWidth',1,'DisplayName','Temperatura superficie')
hold on
%plot(time(jj-10:end),Tfluid(jj-10:end))
plot(time(jj-10:end),Tinterf(jj-10:end),'LineWidth',1,'DisplayName','Temperatura Interfaccia')
grid on
xlabel('tempo [s]')
ylabel('Temperatura [°C]')
legend('Location','best')
title('Evoluzione temperatura nel tempo durante transizione in H-mode')

figure(7)
plot(time(jj-10:end),Ts(jj-10:end),'LineWidth',1,'DisplayName','Temperatura superficie (VS)')
hold on
%plot(time(jj-10:end),Tfluid(jj-10:end))
plot(time(jj-10:end),Tinterf(jj-10:end),'LineWidth',1,'DisplayName','Temperatura Interfaccia (VS)')


figure (5)
subplot(1,2,1)
plot((time(jj+1:end)-t0)*1e3,qs./1e6,'LineWidth',0.5)
grid on
title('Carico dovuto a ELMs')
xlabel('tempo [ms]')
ylabel('q'''' [MW/m^2]')

subplot(1,2,2)
plot(time(jj-10:end)*1e3,Ts(jj-10:end),'LineWidth',1,'DisplayName','Temperatura superficie')
hold on
%plot(time(jj-10:end),Tfluid(jj-10:end))
plot(time(jj-10:end)*1e3,Tinterf(jj-10:end),'LineWidth',1,'DisplayName','Temperatura Interfaccia')
grid on
xlabel('tempo [ms]')
ylabel('Temperatura [°C]')
legend('Location','best')
title('Evoluzione temperatura nel tempo durante H-mode')

figure(10)
plot(time,qev)
grid on
ylabel('\Phi_e_v [W/m2]')
xlabel('Tempo [s]')
%%

T0=T_b4_elms;
dt=.5e-3;
idx=1;
ooo=1;
ppp=1;
istanti=[0.05 0.1 0.151];
for i=1:0.2/dt
    time(jj+i)=time(jj+i-1)+dt;
    if dt*ooo<=deltac
        qs(i)=qmin;
        ooo=ooo+1;
    else
         if dt*ooo-deltac<deltap
            if dt*i<=0.05 || time(jj+i)>=1.2
                qs(i)=qmax1;
            elseif dt*i<=0.1
                qs(i)=qmax2;
            elseif dt*i<0.2
                qs(i)=qmax3;
            end
            ooo=ooo+1;
        else
            qs(i)=qmin;
            ooo=1;
         end
    end
    errc=1;
    
    while errc>toll
        idx=idx+1;
       %assemblaggio matrice
        Tcps=T0(1:n1);
        TBNC=T0(n1+1:end);
        
        %%% cps
        Tcpsp=[Tcps(2:end);Tcps(end)];
        Tcpsm=[Tcps(1);Tcps(1:end-1)];
        kicps=k_CPS(Tcps);
        kipcps=k_CPS(Tcpsp);
        kimcps=k_CPS(Tcpsm);

    
        sub_cps=-dt*1/dx1^2/cv_CPS*(kicps-(kipcps-kimcps)./4);
        main_cps=1+2/dx1^2*dt*kicps./cv_CPS;
        sup_cps=-dt*1/dx1^2/cv_CPS*(kicps+(kipcps-kimcps)./4);
        
        %%% Bulk
        aa=dt*k_BNC/dx2^2/cv_BNC;
        sub_b=-aa*ones(length(TBNC),1);
        main_b=1+2*aa*ones(length(TBNC),1);
        sup_b=sub_b;
        
        %%% dominio intero
        sub=[sub_cps(2:end);sub_b;0];
        main=[main_cps;main_b];
        sup=[0;sup_cps;sup_b(1:end-1)];
    
    
        AA=spdiags([sub,main,sup],-1:1,nn,nn);
        

        qev(jj+i)=0;
        %BC
        AA([1 n1 end],:)=0;
        
        AA(1,1:3)=[-3 4 -1];
        Tp(1)=-(qs(i)-qev(jj+i))*2*dx1/k_CPS(T0(1));
        
        AA(n1,n1-1)=-kicps(end)/k_BNC*dx2/dx1;
        AA(n1,n1)=kicps(end)/k_BNC*dx2/dx1+1;
        AA(n1,n1+1)=-1;
        Tp(n1)=0;

        AA(end,end)=1;
        Tp(end)=330;
    
        T_NoVS=AA\Tp;
        errc=norm(T_NoVS-T0)/norm(T0-T_water);
        T0=T_NoVS;
    end
    if i*dt==istanti(ppp)
        figure(1)
       
        plot(xx*1e3,T_NoVS,'DisplayName',strcat('Profilo dopo ',num2str(i*dt),' s dall''inizio carico pulsato (NO VS)'),LineWidth=1.8)
        hold on
        if ppp<3
            ppp=ppp+1;
        else
            ppp=3;
        end
    end
    Ts(jj+i)=T_NoVS(1);
    Tinterf(jj+i)=T_NoVS(n1);
    Tfluid(jj+i)=T_NoVS(end);
    errt=norm(T_NoVS-Told)/norm(Told-T_water);
    Tp=T_NoVS;
    Told=T_NoVS;
end
figure(1)
plot(xx*1e3,T_NoVS,'DisplayName','Profilo dopo carico pulsato (NO VS)',LineWidth=1.8)
legend()
figure(7)
plot(time(jj+1:end),Ts(jj+1:end),'LineWidth',1,'DisplayName','Temperatura superficie (NO VS)')
hold on
%plot(time(jj-10:end),Tfluid(jj-10:end))
plot(time(jj+1:end),Tinterf(jj+1:end),'LineWidth',1,'DisplayName','Temperatura Interfaccia (NO VS)')
grid on
xlabel('tempo [s]')
ylabel('Temperatura [°C]')
legend('Location','best')
title('Evoluzione temperatura nel tempo durante H-mode')