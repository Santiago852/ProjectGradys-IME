
%% Author: Andrezo
% IME - 08/08/2021
% 
% --------------------V1(08/08/2021)----------------------
% V1 : Two-Ray Model with:
%       * Vertical and Horizontal Polarization
%       * Espherical Earth
%       * Ground Fresnel Reflection Coefficients
%       * Divergence Factor
%       * Ground Spreading Factor
% --------------------V2(11/08/2021)----------------------
% V2 : Two-Ray Model with:
%       * Antennas Profile
%       * Contour Plot of Power Margin
%       * Note 1: Atmospheric attenuation for SHF or higher frequencies
%           (> 10 GHz) not included
%       * Note 2: Diffraction not included
%       * Note 3: Vectorial Sum for Vertical Polarization not included
% --------------------V3(18/08/2021)----------------------
% V3 : Two-Ray Model with:
%       * Polarization Loss Factor (PLF) (use only the same
%           polarization in TX and RX
%       * Incluir Power Delay Profile
%       * Incluir a modelagem do Terreno
%       * Incluir reflexões laterias
% --------------------V4(11/11/2021)----------------------
% V4 : Using some values measured with ESP32 - WiFi, such as:
%       * Real Dipole Antenna
%       * Ground with different soils
%       * Some measured results included
%       * Maximal Two-Ray Model variations included
%       * Note1: For validation of ground parameter, look for:
%           Calc_CoefRealTest7 and forward
% --------------------V3(13/01/2022)----------------------
% V5 : Some propagation errors corrected.

%% Initiation:

close all;
clear all;
format long e

%% Global and Constants:

global Req c f beta sigma;
R = 6371e3;         % Earth Radius (m)
c = 3e8;            % Light Speed (m/s)

%% Data Input:

% Parâmetros Rádio
f = 2.412e3;          % Operational Frequency (MHz).
% Wi-Fi and Bluetooth: f = 2.4e3, usually. Please especify the channel
% frequency. Wi-Fi channels: 1: 2.412, 2: 2.417, 3: 2.422, and so forth
% always increasing in 0.005 GHz for each channel.
beta=2*pi*(f*1e6)/c;% Phase Constant
Pt_dBm = 10;      % TX Power (dBm). ESP32 has some predefined PTXs.
% Put here the chosen value.
Sr_dBm = -85;       % RX Sensitivity. Put here the minimum observed RSSI value.
Gt = 1.97;             % TX Antenna Maximun Gain
Gr = Gt;            % RX Antenna Maximun Gain
Lad1 = 12.3;         % Additional Losses, such as in TX and RX connectors.
% Perdas internas de cada ESP32: 1.5dB. Perda interna do AdalmPluto. 0.3dB.
% Assim, se forem 2 ESP, Lad=3. Se for 1 ESP e 1 Pluto, Lad = 1.8.
% Perda quando ligo o ESP direto no Satsagen e meço o CP do CW: 2.5 dB.
% Perda quando ligo o ESP no Satsagen usando o divisor de potência e meço
% o pico (que deveria ser 10 dBm): entre 11.18 e 12 dB. Na dúvida, colocar 12 dB.
% (Daria um pouco mais de perda, mas considero que o cabo do AdalmPluto deu
% uns 0.3 dB de perda adicional.
Lad2 = 14.22;       % Estava sendo 14.22 até o dia 14/02/22.
% 2º termo (Lad2) observado em medidas. Relativo à mudança de modulação
% (Modulation and Coding Scheme - MSC), mas eu tenho que confirmar isso.
% Put here a value so that Pt_dBm + Gt + Gr - Lad = Maximum Observed Power
pol = 2;            % TX and RX Antenna Polarization (1-Horizontal, 2-Vertical)
ant = 1;            % TX and RX Antenna Type (1-Theoretical Half-Wave Dipole)
                    % (2-Real ESP32 Dipole)

% Parâmetros da terra
K = 1.33;           % Effective Earth Radius Factor:
% K = 1.33 = (4/3) for Standard Atmosphere;  K = inf for plane Earth
Req = K*R;          % Effective Earth Radius
er1 = 1;            % Relative Permittivity - Mean 1 (Air: 1)
er2_1 = 1.7;           % Relative Permittivity - Mean 2 (Poor Ground (dry): 4;
% Average Grounds: 15; Good Ground (wet): 25; Water: 81)
er2_2 = 42;
er_1 = er2_1/er1;       % Relation between relative permittivities 
er_2 = er2_2/er1;
sig = 0;            % Conductivity (mS)(Poor Ground: 1; Average Grouns: 5;
% Good Ground (wet): 20; Destilated Water: 10; Sea Water: 5e3)
sigma = sig*10^(-3);
delta_h_1 = 0.1;     % Ground height variations (cm)
delta_h_1 = delta_h_1*1e-2;
delta_h_2 = 14.7;     % Ground height variations (cm)
delta_h_2 = delta_h_2*1e-2;
ttype = 3;          % 1 for all covered with (er_1,delta_h_1); 2 for all covered with (er_2,delta_h_2)
                    % and 3 for Aterro do Flamengo.

% Parâmetros de posicionamento
ht = 5;            % Transmitter height (m)
hr = 0.15;             % Receiver height (m)
dx = 0:0.1:200;       % Distance between transceptors (m)
dy = 0:0.1:200;       % Distance between transceptors (m)

%% Two-Ray Model Calculation

[d1,d2,ht,hr,r_d,r_r,D,theta,theta_d_tx,theta_d_rx,theta_r,phi]=espacial3(ht,hr,dx,dy);

[G_d,G_r] = gains1(theta_d_tx,theta_d_rx,theta_r,phi,pol,ant,Gt);

[er,delta_h] = terrain(d2,er_1,er_2,delta_h_1,delta_h_2,ttype);

[TwoRL,R1,R2] = two_ray5(theta,pol,delta_h,D,r_d,r_r,G_d,G_r,er);

TwoRLmax = two_ray_simp2(theta,pol,r_d,r_r,G_d,G_r);

% Perda no Espaço Livre (dB). Lembrar que a distância é a composta nos 3
% eixos e deve estar em km
FSL = 20*log10(r_d./1e3)+20*log10(f)+32.45;
% Ajustes de FSL para d=0. Se d=0, calcula p/ a metade do passo.
for i=1:length(dx)
    for k=1:length(dy)
        if (r_d(k,i)==0)
            r_d1 = (r_d(k+1,i+1)/2)/1e3;
            FSL(k,i) = 20*log10(r_d1)+20*log10(f)+32.45;
        end
    end
end

% Att = FSL+TwoRL;                % Atenuação Total
Pr_2R = Pt_dBm-Lad1-Lad2-TwoRL;        % Potência Recebida: Modelo de 2 Raios

%Pr_2R e um vetor 

Pr_2Rmax = Pt_dBm-Lad1-Lad2-TwoRLmax;   % Potência Recebida: Modelo de 2 Raios
% sem levar em consideração os parâmetros do solo, ou seja, raio totalmente
% refletido com inversão de fase. (curva azul ciano)
Pr_FSL = Pt_dBm+G_d-Lad1-Lad2-FSL;      % Potência Recebida: Espaço Livre 
Pr_ISO = Pt_dBm+Gt+Gr-Lad1-Lad2-FSL;    % Potência Recebida: Isotrópica. Sem levar
% consideração os diagramas de irradiação.
Pr_R1 = Pt_dBm-Lad1-Lad2+20*log10(abs(R1(:,:,1)));
Pr_R2 = Pt_dBm-Lad1-Lad2+20*log10(abs(R2(:,:,1)));

% Margens de Cobertura
M_2R = Pr_2R - Sr_dBm;
M_FSL = Pr_FSL - Sr_dBm;
M_ISO = Pr_ISO - Sr_dBm;
M_R1 = Pr_R1 - Sr_dBm;
M_R2 = Pr_R2 - Sr_dBm;

% Ajuste p/ não mostrar potência onde não vai funcionar. Também calcula a
% área de cobertura.
C_2R = 0;
C_FSL = 0;
C_ISO = 0;
C_R1 = 0;
C_R2 = 0;
dy_max_2R = 0;
dy_max_FSL = 0;
dy_max_ISO = 0;
PDP1 = [];
Delay1 = [];
PDP2 = [];
Delay2 = [];

for i=1:length(dx)
    for k=1:length(dy)
        if (M_2R(k,i)<0)
            M_2R(k,i) = -inf;            
        else
            C_2R = C_2R+1;  % Contagem de pontos com cobertura.
            
            % Cálculo da distância máxima de cobertura
            if (dy(k)>dy_max_2R)
                dy_max_2R = dy(k);
            end
        end
        if (M_FSL(k,i)<0)
            M_FSL(k,i) = -inf;            
        else
            C_FSL = C_FSL+1;  % Contagem de pontos com cobertura.
            
            % Cálculo da distância máxima de cobertura
            if (dy(k)>dy_max_FSL)
                dy_max_FSL = dy(k);
            end
        end
        if (M_ISO(k,i)<0)
            M_ISO(k,i) = -inf;
        else
            C_ISO = C_ISO+1;  % Contagem de pontos com cobertura.
            
            % Cálculo da distância máxima de cobertura
            if (dy(k)>dy_max_ISO)
                dy_max_ISO = dy(k);
            end
        end
        if (M_R1(k,i)>= 0)
            % Power Delay Profile (Raio Direto)
            C_R1 = C_R1+1;  % Contagem de pontos p/ o PDP do Raio Direto.            
            PDP1(C_R1) = M_R1(k,i);
            Delay1(C_R1) = R1(k,i,3);
        end
        if (M_R2(k,i)>= 0)
            % Power Delay Profile (Raio Refletido)
            C_R2 = C_R2+1;  % Contagem de pontos p/ o PDP do Raio Refletido.            
            PDP2(C_R2) = M_R2(k,i);
            Delay2(C_R2) = R2(k,i,3);
        end
    end
end

[Delay1,PDPSort] = sort(Delay1);
PDP1 = PDP1(PDPSort);
[Delay2,PDPSort] = sort(Delay2);
PDP2 = PDP2(PDPSort);

TotArea = (max(dx)-min(dx))*(max(dy)-min(dy));
Npoints = length(dx)*length(dy);
% Cov_2R = 4*C_2R*TotArea/Npoints;
% Cov_FSL = 4*C_FSL*TotArea/Npoints;
% Cov_ISO = 4*C_ISO*TotArea/Npoints;
Cov_2R = C_2R*100/Npoints;
Cov_FSL = C_FSL*100/Npoints;
Cov_ISO = C_ISO*100/Npoints;

%% Measured Results

% Resultdos de Testes de Propagação feitos no Aterro em 10/12/21, com a
% Amanda e consolidados por mim.
% Cada linha uma distância diferente, mas sempre com hr=0.15m e ht=5m.

% Resultados medidos no lóbulo principal (theta=0º) da antena dipolo na
% horizontal.
dist = [0 20 40 50 60 80 83];
Mean0 = [-57.7 -85.725 -85.85 -88.75 -95.775 -94.4 -97.8];
Max0 = [-57.2 -82.2 -84.7 -88.6 -94.4 -92.9 -92.8];

% Resultdos de Testes de Propagação feitos no Aterro em 15/12/21, com o
% Pedro e consolidados por mim. Peguei as duas meédias dos arquivos e fiz a
% média dividindo a soma por 2. Incluí também os resultados p/ d = 0; 40 e
% 50, feitos com a Amanda no dia anterior.
% Incluí também os resultados das medições no dia 21/12/21 p/ as distâncias
% de 0.5 1 e 2m. Daí preferi tirar o resultado anterior p/ 0m.
% Cada linha uma distância diferente, mas sempre com hr=0.15m e ht=5m.

% Resultados medidos no lóbulo principal (theta=0º) da antena dipono na
% horizontal.
dist1 = [0.5 1 2 2.2 3.5 6 11.5 18 20 30 40 50];
Mean1 = [-62.95 -66.7 -72.35 -75.4 -71.4 -70.05 -80.8 -76.3 -77.95 -84.55 -85.85 -88.75];

% Resultados medidos no lóbulo principal (theta=0º) da antena dipono na
% vertical.
dist2 = [0.4 2.2 3.5 6 10.5 20 27.5 40];
% Medidas anteriores, feitas em 28/12/21.
% Mean02 = [-79.55 -74.5 -67.5 -77.55 -71.95 -75.25 -79.75 -81.65]; %
Mean2 = [-79.55 -82.1 -72.75 -77.55 -71.95 -75.25 -79.75 -81.65];

dist5= [2.2 6 10.5 20 27.5 40];
Mean5_semAtenuacao= [-65.9339 -58.8057 -58.2252 -59.5710 -61.0619 -62.6337];
atenuacao5=[-0.52375 -0.449 -0.49 -0.5675 -0.788888889 -1.1585];
Mean5=Mean5_semAtenuacao+atenuacao5

% Medidas feitas com sinal CW em 14/02/22 pegando o MaxHold e uma média
% entre o máximo e o mínimo.
dist3 = [2 3.5 6 11.5 18 30];
Max1 = [-50 -49 -49 -55.5 -53 -60];
Mean3 = [-53 -50 -50 -62 -56 -67];

% Medidas feitas com sinal CW em 05/05/22 pegando o pico do Satsagen em 2.4
dist4 = [2.0000    3.5000    6.0000   11.5000   18.0000   30.0000   40.0000];
Mean4_semAtenuacao = [-53.2916  -62.7664  -53.0484  -66.6641  -55.9641  -58.3741  -60.7583];
atenuacao4=[-0.053636364 -0.097272727 -0.046363636 -0.01375 0.09875 0.08 -0.54];
Mean4 = Mean4_semAtenuacao+atenuacao4;

%%COLOCAR TAMBÉM AS MEDIDAS DE CP DO SATSAGEN!!!!!
 
%% Plot Image

figure
hold on;
plot(dy,Pr_2R(:,1),'b','Linewidth',2);
% b = [Pr_2R(:,1)-min(Pr_2Rmax(:,1),Pr_FSL(:,1)) max(Pr_2Rmax(:,1),Pr_FSL(:,1))-Pr_2R(:,1)];
% [l,p] = boundedline(dy,Pr_2R(:,1),b,'-','alpha');
% outlinebounds(l,p);

%--plot(dy,Pr_2Rmax(:,1),'c','Linewidth',2);
%--plot(dy,Pr_FSL(:,1),'--b','Linewidth',2);
%--plot(dy,Pr_2R(:,1)+Lad2,'k',dy,Pr_FSL(:,1)+Lad2,'--k','Linewidth',2);

% plot(dist,Mean0,'-ob','Linewidth',2);



if pol == 1
    %plot(dist1,Mean1,'-or','Linewidth',2);
    % plot(dist,Max0,'-og','Linewidth',2);
    %plot(dist3,Mean3,'-og','Linewidth',2);
    
    %--plot(dist4,Mean4,'-oy','Linewidth',2);
    
    %pegar os pontos da curva preta nos pontos da dist do experimento do Aterro
%curva preta Pr_2R(:,1)+Lad2
dy_y=reshape(dist4,[1,7]);
for i=1:7
    for r=1:501
    if dy_y(i)==dy(r)
       Pr_2R(r,1)+Lad2
Delta(i)=Mean4(i)-(Pr_2R(r,1)+Lad2)

    end
    end
end
%--plot(dy_y,Delta,'m')
elseif pol == 2
    
    %plot(dist2,Mean2,'-or','Linewidth',2);
    %--plot(dist5,Mean5,'-oy','Linewidth',2);
 dx_x=reshape(dist5,[1,6]);
    for i=1:6
    for r=1:501
    if dx_x(i)==dx(r)
       Pr_2R(r,1)+Lad2
Delta(i)=Mean5(i)-(Pr_2R(r,1)+Lad2)

    end
    end
end
%--plot(dx_x,Delta,'m')
end

%--plot(dy,Sr_dBm*ones(length(dy)),'k','Linewidth',2);
xlabel('dx (m)');
ylabel('Pr(dBm)');
grid
hold off;



% figure
% b = [Pr_2R(:,1)-min(Pr_2Rmax(:,1),Pr_FSL(:,1)) max(Pr_2Rmax(:,1),Pr_FSL(:,1))-Pr_2R(:,1)];
% [l,p] = boundedline(dy,Pr_2R(:,1),b,'-','alpha');
% outlinebounds(l,p);
% xlabel('dx (m)');
% ylabel('Pr(dBm)');
% grid

% figure
% contourf(dx,dy,M_2R);
% title('Two Ray Model Coverage');
% xlabel('dx (m)');
% ylabel('dy (m)');
% zlabel('M(dB)');
% colorbar
% legend('Margin above sensitivity (dB)');
% grid
% dim = [0.48 0.73 0.1 0.1];
% str = sprintf('TX Height: %4.1f m\nRX Height: %4.1f m\nCoverage: %4.2f %%\nMax Distance: %4.1f m', ht, hr,Cov_2R,dy_max_2R);
% annotation('textbox',dim,'String',str,'FitBoxToText','on');
% % text(350,400,str);
% 
% figure
% contourf(dx,dy,M_FSL);
% title('Free Space Loss Model Coverage');
% xlabel('dx (m)');
% ylabel('dy (m)');
% zlabel('M(dB)');
% colorbar
% legend('Margin above sensitivity (dB)');
% grid
% str = sprintf('TX Height: %4.1f m\nRX Height: %4.1f m\nCoverage: %4.2f %%\nMax Distance: %4.1f m', ht, hr,Cov_FSL,dy_max_FSL);
% annotation('textbox',dim,'String',str,'FitBoxToText','on');

% figure
% contourf(dx,dy,M_ISO);
% title('Isotropic Coverage');
% xlabel('dx (m)');
% ylabel('dy (m)');
% zlabel('M(dB)');
% colorbar
% legend('Margin above sensitivity (dB)');
% grid
% str = sprintf('TX Height: %4.1f m\nRX Height: %4.1f m\nCoverage: %4.2f %%\nMax Distance: %4.1f m', ht, hr,Cov_ISO,dy_max_ISO);
% annotation('textbox',dim,'String',str,'FitBoxToText','on');

% figure
% hold on;
% plot(Delay1,PDP1,Delay2,PDP2,'Linewidth',2);
% title('Power Delay Profile');
% xlabel('Delay (s)');
% ylabel('Received Power');
% grid
% hold off;

% Figuras p/ mostrar as permissividades e delta_hs considerados para cata
% TX na superfície.
% figure
% h = surf(dx,dy,er);
% set(h,'LineStyle','none')
% 
% figure
% h = surf(dx,dy,delta_h);
% set(h,'LineStyle','none')
