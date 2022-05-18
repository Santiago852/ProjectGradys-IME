function [TwoRL] = two_ray_simp2(theta,pol,r_d,r_r,G_d,G_r)
%% Calcula a atenuação pelo modelo de 2 raios
% ht e hr em metros, d em km
% Considera o raio totalmente refletido com inversão de fase, ou seja,
% Gamma = -1, Rho_s = 1 e D = 1.
% Author: Andrezo
% V1: 11/08/2021

%% Initiation:
    global beta c;
            
%% Calculation of Attenuation
    % Coeficiente de Reflexão de Fresnel: F = [Magnitude Fase(rad)]
    % F = fresnel2(theta,pol);
    
    % Fator de Perda de Espalhamento
    % delta_fase = 2*beta*delta_h.*cos(theta);
    % rho_s = exp(-0.5.*(delta_fase.^2));
    
    % Conversão dos Ganhos dos Raios Direto e Refletido
    F_d = sqrt(10.^(G_d./10));
    F_r = sqrt(10.^(G_r./10));
    
    % Atenuação pelo Modelo de 2 Raios (dB)
    
    % Cálculo do Raio Direto já com a FSL só dele
    r1 = F_d.*exp(-j*beta.*r_d)./(2*beta.*r_d);
        
    % Cálculo do Raio Refletido já com a FSL só dele
    r2 = (-1)*F_r.*exp(-j*beta.*r_r)./(2*beta.*r_r);
        
    L = abs(r1+r2);
    TwoRL = -20*log10(L);
end