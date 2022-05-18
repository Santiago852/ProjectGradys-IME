function [G_d,G_r] = gains1(theta_d_tx,theta_d_rx,theta_r,phi,pol,ant,Gmax)
%% Calcula os ganhos direcionais de cada antena no sistema
% l é o comprimento elétrico (em comprimentos de onda); theta em radianos;
% Gmax em dBi e pol ainda está como opção (1 - vertical, 2 - horizontal).
% Os resultados de G_d e G_r são os ganhos de campo em relação direta.
% Author: Andrezo
% V1: 12/08/2021

%% Calculation of Gains  
switch pol
    case 1
        disp('A polarizaçao da onda refletida no solo é: Horizontal.');
        % Nessa configuração o TX está no drone com o seu plano
        % omnidirecional sempre apontado na direção do RX.
        % O RX está na origem com a antena alinhada com x (plano
        % ominidirecional em y-z, logo o theta_tx será pi/2 e o theta_rx
        % é o acos(sin(theta).*cos(phi), o que é meio complicado de ver,
        % mas confirmei escrevendo isso.
        theta_d_tx = pi/2;
        theta_r_tx = pi/2;
        theta_d_rx = acos(sin(theta_d_rx).*cos(phi));
        theta_r_rx = acos(sin(theta_r).*cos(phi));
        % Inclusão do PLF para a pol. horizontal.
        PLF = sin(phi).^2;
        PLF = 10*log10(PLF);
    case 2
        disp('A polarizaçao da onda refletida no solo é: Vertical.');
        theta_r_tx = theta_r;
        theta_r_rx = theta_r;
        % Inclusão do PLF para a pol. vertical. É sempre 1.
        PLF = 1;
        PLF = 10*log10(PLF);
    otherwise
        disp('Opção inválida');
    return
end

switch ant
    case 1
        disp('A antena é uma dipolo de meia onda.');
        l = 1/2;
        
        G_d_tx = vdipole1(l,theta_d_tx,Gmax);
        G_d_rx = vdipole1(l,theta_d_rx,Gmax);
        G_d = G_d_tx+G_d_rx+PLF;
                
        G_r_tx = vdipole1(l,theta_r_tx,Gmax);
        G_r_rx = vdipole1(l,theta_r_rx,Gmax);
        G_r = G_r_tx+G_r_rx+PLF;
    
    case 2
        disp('A antena é a dipolo de meia onda do ESP32.');
        disp('Obs: O TX, por ser considerado acima do RX, teve a sua antena virada pra baixo.');        
        G_d_tx = realvdipole(theta_d_tx,Gmax);
        G_d_rx = realvdipole(theta_d_rx,Gmax);
        G_d = G_d_tx+G_d_rx+PLF;
                
        G_r_tx = realvdipole(theta_r_tx,Gmax);
        G_r_rx = realvdipole(theta_r_rx,Gmax);
        G_r = G_r_tx+G_r_rx+PLF;
        
    otherwise
        disp('Opção inválida');
    return
end

end