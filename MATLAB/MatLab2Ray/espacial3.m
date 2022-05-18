function [d1,d2,ht,hr,r_d,r_r,D,theta,theta_d_tx,theta_d_rx,theta_r,phi]  = espacial3(ht_,hr_,dx_,dy_)
%% Calcula os par�metros espaciais para o modelo de 2 raios
% ht, hr, dx e dy em metros
% Author: Andrezo
% V1: 11/08/2021

%% Initiation:
    global Req;
    
%% Calculation of Spatial Parameters
    % Criando a dist�ncia equivalente em cada dx e dy
    dx = repmat(dx_,length(dy_),1);
    dy = repmat(dy_',1,length(dx_));
    d = sqrt(dx.^2+dy.^2);  % Dist�ncia equivalente em cada dx e dy            
    dh = sqrt(2*Req*max(ht_,hr_));    % Maior dist�ncia p/ o Horizonte-R�dio
    if (d<=dh)                      % Teste de Terra Plana
        fprintf('Terra Plana. Dist�ncia para o horizonte: %5.2f km.\n', dh/1000);
        
        %Dist�ncias p/ o ponto de reflex�o na superf�cie (TX e RX) (m)
        d1 = ht_.*d./(ht_+hr_);
        d2 = hr_.*d./(ht_+hr_);
        
        % Alturas Equivalentes do TX e do RX (m)
        ht = ht_;
        hr = hr_;
        
        % Trajeto do raio direto e diferen�a do refletido p/ o direto (m)
        r_d = sqrt(d.^2+(ht-hr).^2); 
        r_r = sqrt(d.^2+(ht+hr).^2);
        % dr = 2*ht*hr./d;
        dr = abs(r_r-r_d);  % Melhor usar essa defini��o pra evitar um
        % infinito que aparecia quando d=0.
        
        % Fator de Diverg�ncia
        D = ones(size(d));
        
    else
        fprintf('Terra Esf�rica. Dist�ncia para o horizonte: %5.2f km.\n', dh/1000);
        
        %Dist�ncias p/ o ponto de reflex�o na superf�cie (TX e RX) (m)
        p = 2*sqrt((Req*(ht_+hr_)+d.^2./4)./3);
        phi = acos(2*Req*abs(ht_-hr_)*d./p.^3);
        d1 = d./2+p.*cos((phi+pi)/3);
        d2=hr.*d1/ht;
    
        % Alturas Equivalentes do TX e do RX (m)
        ht = ht_-(d1.^2)./(2*Req);
        hr = resolver((d1./ht).^2,2*Req,-2*Req*hr_); % a fun��o resolver sempre 
        % vai retornar o valor positivo para a equa��o do 2� grau.
    
        % Trajeto do raio direto e diferen�a do refletido p/ o direto (m)
        r_d = sqrt(d.^2+(ht-hr).^2); 
        r_r = sqrt(d.^2+(ht+hr).^2);
        % dr = 2*ht*hr./d;
        dr = abs(r_r-r_d);  % Melhor usar essa defini��o pra evitar um
        % infinito que aparecia quando d=0.
    
        % Fator de Diverg�ncia
        D = (1+(2.*d1.*d2)./(Req*(ht+hr))).^(-1/2);
        
    end
    
    % �ngulo de Incid�ncia no solo
    theta = atan(d1./ht);
    
    % �ngulos de Azimute (em rela��o a z) do raio direto (TX e RX)
    % Como atan � negativo p/ tangentes negativas, podemos fazer uma s�
    % express�o mesmo que ht<hr. �ngulos em rela��o a z positivo.
    theta_d_tx = pi/2+atan((ht-hr)./d);
    theta_d_rx = pi/2-atan((ht-hr)./d);
    
    % �ngulos de Azimute (em rela��o a z) do raio direto (TX e RX)
    theta_r = pi/2+atan((ht+hr)./d);
    % Como os thetas pro raio refletido em TX e RX s�o iguais, usei s� uma
    % vari�vel, ent�o. � bom ter o desenho dos raios direto e refletido pra
    % entender essas express�es.
    
    % �ngulo de Azimute (em rela��o a x) dos raios direto e refletido
    % (RX colocado no centro).
    phi = atan(dy./dx);
    
    % Ajustes dos �ngulos para d=0.
    for i=1:length(dx_)
        for k=1:length(dy_)
            if (d(k,i)==0)
                phi(k,i) = pi/2; % Vulgo 45�.
                if (ht==hr)
                    theta_d_tx(k,i) = pi/2; % Se estiverem na mesma altura,
                    theta_d_rx(k,i) = pi/2; % theta_d = 90�.
                    theta_r(k,i) = pi;
                end
            end
        end
    end

end