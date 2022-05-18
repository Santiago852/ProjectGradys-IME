function [TwoRL] = two_ray_simp2(theta,pol,r_d,r_r,G_d,G_r)
%% Calcula a atenua��o pelo modelo de 2 raios
% ht e hr em metros, d em km
% Considera o raio totalmente refletido com invers�o de fase, ou seja,
% Gamma = -1, Rho_s = 1 e D = 1.
% Author: Andrezo
% V1: 11/08/2021

%% Initiation:
    global beta c;
            
%% Calculation of Attenuation
    % Coeficiente de Reflex�o de Fresnel: F = [Magnitude Fase(rad)]
    % F = fresnel2(theta,pol);
    
    % Fator de Perda de Espalhamento
    % delta_fase = 2*beta*delta_h.*cos(theta);
    % rho_s = exp(-0.5.*(delta_fase.^2));
    
    % Convers�o dos Ganhos dos Raios Direto e Refletido
    F_d = sqrt(10.^(G_d./10));
    F_r = sqrt(10.^(G_r./10));
    
    % Atenua��o pelo Modelo de 2 Raios (dB)
    
    % C�lculo do Raio Direto j� com a FSL s� dele
    r1 = F_d.*exp(-j*beta.*r_d)./(2*beta.*r_d);
        
    % C�lculo do Raio Refletido j� com a FSL s� dele
    r2 = (-1)*F_r.*exp(-j*beta.*r_r)./(2*beta.*r_r);
        
    L = abs(r1+r2);
    TwoRL = -20*log10(L);
end