function Rf = fresnel3(theta,pol,er)
%Esta função pega o modulo e o angulo do coeficiente de reflexao de fresnel
%Rf(1)=modulo e Rf(2)=fase(em radianos)
    global c f sigma;
    erc=er-j*sigma*60*c/f;
    
    switch pol
        case 1
            F=(cos(theta)-sqrt(erc).*sqrt(1-(1./erc).*sin(theta).^2))./(cos(theta)+sqrt(erc).*sqrt((1-(1./erc).*sin(theta).^2)));
            Rf(:,:,1)=abs(F);
            Rf(:,:,2)=angle(F);
        case 2
            F=(-cos(theta)+sqrt(1./erc).*sqrt(1-(1./erc).*sin(theta).^2))./(cos(theta)+sqrt((1./erc)).*sqrt(1-(1./erc).*sin(theta).^2));
            Rf(:,:,1)=abs(F);
            Rf(:,:,2)=angle(F);
        otherwise
            disp('opção invalida');
            return
    end
end