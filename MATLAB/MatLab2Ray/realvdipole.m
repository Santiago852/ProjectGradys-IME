function Gain = realvdipole(theta_i,Gmax)
%% Calcula o ganho do dipolo da ESP32 em determinada direção
% theta_i em radianos; Gmax em dBi.
% pol ainda está como opção (1 - vertical, 2 - horizontal)
% Author: Andrezo
% V1: 06/12/2021

%% Calculation of Gain  
    filename = 'Ensaio de Antena-IME-EdicaoAndrezo5-211205.xlsx';
    sheet = 'Cálculos-Diagrama de radiação-P';
    xlRangeXY = 'D5:D365';
    xlRangeYZ = 'F5:F365';

    % Valores obtigos para os ganhos da antena dipolo real do ESP32, desde
    % phi = 0 a 360 e theta = 0 a 360 com variação de 1 grau.
    thetad = 0:1:90;
    DipoloXY = xlsread(filename,sheet,xlRangeXY); % Diagrama Azimutal
    DipoloYZ = xlsread(filename,sheet,xlRangeYZ); % Diagrama de Elevação
    
    % Como só medidos theta de 0 a 90 graus, vamos pegar somente o 1º
    % quadrante
    YZUp(1,:) = DipoloYZ(1:91)';
    YZUp(2,:) = DipoloYZ(361:-1:271)';
    DipoloYZUp = mean(YZUp);

    thetai = theta_i*180/pi; % Conversão para graus
    % Como essa antena tem variações entre o 1º e o 2º quadrantes, fizemos
    % como se ela estivesse virada pra baixo, se theta for maior que 90º,
    % ou seja, puxamos pro 1º quadante.
    if (thetai > 90)
        thetai = 180 - thetai;
    end
    % Se o valor de theta não estiver dentro de 1º de variação:
    Gain = interp1(thetad,DipoloYZUp,thetai) + Gmax;
    % Acha o valor interpolado de DipooloYZ, baseado no valor de thetai.
    
end