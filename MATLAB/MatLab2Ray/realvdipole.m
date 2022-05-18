function Gain = realvdipole(theta_i,Gmax)
%% Calcula o ganho do dipolo da ESP32 em determinada dire��o
% theta_i em radianos; Gmax em dBi.
% pol ainda est� como op��o (1 - vertical, 2 - horizontal)
% Author: Andrezo
% V1: 06/12/2021

%% Calculation of Gain  
    filename = 'Ensaio de Antena-IME-EdicaoAndrezo5-211205.xlsx';
    sheet = 'C�lculos-Diagrama de radia��o-P';
    xlRangeXY = 'D5:D365';
    xlRangeYZ = 'F5:F365';

    % Valores obtigos para os ganhos da antena dipolo real do ESP32, desde
    % phi = 0 a 360 e theta = 0 a 360 com varia��o de 1 grau.
    thetad = 0:1:90;
    DipoloXY = xlsread(filename,sheet,xlRangeXY); % Diagrama Azimutal
    DipoloYZ = xlsread(filename,sheet,xlRangeYZ); % Diagrama de Eleva��o
    
    % Como s� medidos theta de 0 a 90 graus, vamos pegar somente o 1�
    % quadrante
    YZUp(1,:) = DipoloYZ(1:91)';
    YZUp(2,:) = DipoloYZ(361:-1:271)';
    DipoloYZUp = mean(YZUp);

    thetai = theta_i*180/pi; % Convers�o para graus
    % Como essa antena tem varia��es entre o 1� e o 2� quadrantes, fizemos
    % como se ela estivesse virada pra baixo, se theta for maior que 90�,
    % ou seja, puxamos pro 1� quadante.
    if (thetai > 90)
        thetai = 180 - thetai;
    end
    % Se o valor de theta n�o estiver dentro de 1� de varia��o:
    Gain = interp1(thetad,DipoloYZUp,thetai) + Gmax;
    % Acha o valor interpolado de DipooloYZ, baseado no valor de thetai.
    
end