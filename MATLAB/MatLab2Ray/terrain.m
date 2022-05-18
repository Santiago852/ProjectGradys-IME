function [er,delta_h] = terrain(d2,er_1,er_2,delta_h_1,delta_h_2,ttype)
%% Preenche a permissividade relativa e as diferenças de altura do solo de acordo com as posições medidas no Aterro
% Author: Andrezo
% V1: 07/12/2021

%% Preenchimendo das matrizes
    % Criação das matrizes er e delta_h com os valores default do cimento
    switch ttype
        case 1
            disp('A superfície é o cimento.');
            er = ones(size(d2))*er_1;
            delta_h = ones(size(d2))*delta_h_1;
        case 2
            disp('A superfície é a grama.');
            er = ones(size(d2))*er_2;
            delta_h = ones(size(d2))*delta_h_2;
        case 3
            disp('A superfície é o campo do Aterro.');
            er = ones(size(d2))*er_1;
            delta_h = ones(size(d2))*delta_h_1;
            % Mudança dos valores onde temos grama. Reparar que a geometria não é
            % exatamente a do Aterro. Estou fazendo só círculos ao redor do centro
            % mais próximo do monumento.
            I = find(d2>5.1 & d2<=15.6);
            er(I) = er_2;
            delta_h(I) = delta_h_2;
    
            I = find(d2>20.7 & d2<=24.5);
            er(I) = er_2;
            delta_h(I) = delta_h_2;
    
            I = find(d2>28.8 & d2<=44.9);
            er(I) = er_2;
            delta_h(I) = delta_h_2;
    
            I = find(d2>55.1 & d2<=71.2);
            er(I) = er_2;
            delta_h(I) = delta_h_2;
    
            I = find(d2>75.4 & d2<=80.4);
            er(I) = er_2;
            delta_h(I) = delta_h_2;
        otherwise
            disp('opção invalida');
            return
        end
    
end