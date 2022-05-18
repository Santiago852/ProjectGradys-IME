function Gain = vdipole1(l,theta_i,Gmax)
%% Calcula o ganho do dipolo em determinada direção
% l é o comprimento elétrico (em comprimentos de onda); theta em radianos;
% Gmax em dBi e pol ainda está como opção (1 - vertical, 2 - horizontal)
% Author: Andrezo
% V1: 12/08/2021

%% Calculation of Gain  
    x = pi*l;

    Un = ((cos(x.*cos(theta_i)) - cos(x))./(sin(theta_i).*(1-cos(x)))).^2;
    % Eq 4.64 do Balanis já normalizada pelo valor máximo (1-cos(x)^2, como nos
    % slides da minha aula
    % Un = Un./max(Un); % Normalizando Un, pois a anterior não dá muita precisão.
    Gain = 10.*log10(Un) + Gmax;
    % Ganho da dipolo de meia-onda Gmax ~ 2.13 dBi (linha de 50 ohms de cobre).

end