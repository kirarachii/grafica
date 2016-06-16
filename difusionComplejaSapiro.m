function I = difusionColorSapiro(I, dt, numIt, k)
    [xDim, yDim, numCh] = size(I);
    ep = 0.01;
    %indices
    x1 = [1 1:xDim-1];
    x2 = [2:xDim xDim];
    y1 = [1 1:yDim-1];
    y2 = [2:yDim yDim];
    for it=1: numIt
        %derivada en direccion x
        Ix1 = I - I(x1,:,:);
        Ix2 = I(x2,:,:) - I;
        Ix = (Ix1 + Ix2)./2;
        %derivada en direccion y
        Iy1 = I - I(:,y1,:);
        Iy2 = I(:,y2,:) - I;
        Iy = (Iy1 + Iy2)./2;
        %sumamos todas las componentes por canal de tensor de estructura
        g11 = sum(Ix.*Ix,3);
        g12 = sum(Ix.*Iy,3);
        g22 = sum(Iy.*Iy,3);
        %calculamos eiganvalues
        l1 = 0.5.*(g11+g22 + sqrt((g11-g22).^2 +4.*g12.^2));
        l2 = 0.5.*(g11+g22 - sqrt((g11-g22).^2 +4.*g12.^2));
        %angulo minima y maxima variacion
        tetap = 0.5.*(atan((2.*g12)./(g11-g22+ep)));
        tetap2 = tetap + pi;
        tetam = tetap + pi/2;
        %componentes x e y de la nueva base vectorial
        etax = cos(tetap);
        etay = sin(tetap);
        xix = cos(tetam);
        xiy = sin(tetam);
        
        %segunda derivada
        Ixx = I(x1,:,:) + I(x2,:,:) -2.*I;
        Iyy = I(:,y1,:) + I(:,y2,:) -2.*I;
        Ixy = 0.25.*(I(x1,y1,:) + I(x2,y2,:) -I(x1,y2,:) - I(x2,y1,:));
        %derivada direccional en direccion a la minima variacion
        Ixixi = xix(:,:,ones(numCh,1)).^2.*Ixx + 2.*xix(:,:,ones(numCh,1)).*xiy(:,:,ones(numCh,1)).*Ixy + xiy(:,:,ones(numCh,1)).^2.*Iyy;
        N = sqrt(l1+l2);
        G = g(N,k);
        I = I + dt.*G(:,:,ones(numCh,1)).*Ixixi;

    end
end

function pondEspacial = g(N,k)
    pondEspacial = 1./(1 + (N./k).^2);
end

%para llamar la funcion Isapiro = difusionColorSapiro(I, 0.1, 100, 5);
