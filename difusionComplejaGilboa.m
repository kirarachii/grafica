function I = difusionComplejaGilboa(I, numIt, dt, k)
    [xDim, yDim] = size(I);
    teta = pi./1000;
    %indices
    xp = [1 1:xDim-1];
    xs = [2:xDim xDim];
    ye = [1 1:yDim-1];
    yo = [2:yDim yDim];
    
    for it = 1:numIt
        I2 = I;
        I = imfilter(I, fspecial('gaussian', [3 3], 1), 'conv', 'symmetric', 'same');
        Ixn = I(xp,:) - I;
        Ixs = I(xs,:) - I;
        Iye = I(:,ye) - I;
        Iyo = I(:,yo) - I;
        
        gn = exp(sqrt(-1).*teta)./(1+(imag(Ixn)./(k.*teta).^2));
        gs = exp(sqrt(-1).*teta)./(1+(imag(Ixs)./(k.*teta).^2));
        ge = exp(sqrt(-1).*teta)./(1+(imag(Iye)./(k.*teta).^2));
        go = exp(sqrt(-1).*teta)./(1+(imag(Iyo)./(k.*teta).^2));
        
        I = I2;
        Ixn = I(xp,:) - I;
        Ixs = I(xs,:) - I;
        Iye = I(:,ye) - I;
        Iyo = I(:,yo) - I;

        I = I + dt.(Ixn.*gn + Ixs.*gs + Iye.*ge + Iyo.*go);
    end
    I = real(I);
end
%para llamarlo Ig = difusionComplejaGilboa(I, 10, 0.1, 5)
%In = imnoise(uint8(I),'gaussian', 0, 0.01);
%figure, imagesc(In), colormap(gray);
%Ig = difusionComplejaGilboa(double(In), 30, 0.1, 5); pueden cambiar las iteracion
%figure, imagesc(uint8(Ig)), colormap(gray);

