function [TwoRL,R1,R2] = two_ray5(theta,pol,delta_h,D,r_d,r_r,G_d,G_r,er)
%% Calcula a atenuação pelo modelo de 2 raios
% ht e hr em metros, d em km
% Author: Andrezo
% V1: 11/08/2021

%% Initiation:
    global beta c;
            
%% Calculation of Attenuation
    % Coeficiente de Reflexão de Fresnel: F = [Magnitude Fase(rad)]
    F = fresnel3(theta,pol,er);
    
    % Fator de Perda de Espalhamento
    delta_fase = 2*beta.*delta_h.*cos(theta);
    rho_s = exp(-0.5.*(delta_fase.^2));
    
    % Conversão dos Ganhos dos Raios Direto e Refletido
    F_d = sqrt(10.^(G_d./10));
    F_r = sqrt(10.^(G_r./10));
    
    % Atenuação pelo Modelo de 2 Raios (dB)
    
    % Cálculo do Raio Direto já com a FSL só dele
    r1 = F_d.*exp(-j*beta.*r_d)./(2*beta.*r_d);
    R1(:,:,1) = abs(r1);            % Amplitude do Rd
    R1(:,:,2) = angle(r1);          % Fase do Rd
    R1(:,:,3) = r_d/c;              % Retardo do Rd
    
    % Cálculo do Raio Refletido já com a FSL só dele
    r2 = F_r.*F(:,:,1).*D.*rho_s.*exp(-j*(beta.*r_r+F(:,:,2)))./(2*beta.*r_r);
    R2(:,:,1) = abs(r2);            % Amplitude do Rr
    R2(:,:,2) = angle(r2);          % Fase do Rr
    R2(:,:,3) = r_r/c;              % Retardo do Rr
    
    L = abs(r1+r2);
    TwoRL = -20*log10(L);
end