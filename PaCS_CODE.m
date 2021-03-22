clc
clear
pois=0.48;%poisonss ratio
young=2e3;%youngs modulus
young_e12= young; %* 1e-12;
pixelsize=.227;

Af=[2144.33];%final area of pattern in Um^2
for kk=1:length(Af)
A_f=Af(kk);
A_i=2401.96; %initial area of pattern in Um^2
%%
bound=sqrt(A_i)/2;
bound2=sqrt(A_f)/2;
%NN=20-1;
%kk=1;
NN=20;%grid size
%for NN=4:2:200
    NN=NN-1;
    int=2*bound/NN;
     X1=-bound:int:bound;
     %X(length(X))=bound;
     Y1=-bound:int:bound;
     [x1,y1]=meshgrid(X1,Y1); 

    int2=2*bound2/NN;
     X2=-bound2:int2:bound2;
     %X(length(X))=bound;
     Y2=-bound2:int2:bound2;
     [x2,y2]=meshgrid(X2,Y2);

     ux=x2-x1;
     uy=y2-y1;

    %%
    a1   = (1.0 + pois) * (1.0 - pois) / (pi * young_e12);
    b1   = (1.0 + pois) * pois / (pi * young_e12);
    c1   = (1.0 + pois) * pois / (pi * young_e12);
    %%
    % x_org = x;
    %     y_org = y;
    %     
    %     x  = (x - mean(x(:))) * pixelsize;
    %     y  = (y - mean(y(:))) * pixelsize;
    % %     ux = dx * pixelsize;
    % %     uy = dy * pixelsize;

        xv = x1(:);
        yv = y1(:);
        xv_orig = xv;
        yv_orig = yv;

        spacing = yv(2) - yv(1);

        N   = sqrt(length(xv));

        %%
        %  kx, ky, AND  k_abs

        clear mala*;
        for i = 1 : (N/2)
            malax(i,:) = 0 : ((N/2)-1);
            malay(i,1:(N/2)) = (N/2)-i;
        end

        kx = [ malax  malax-(N/2);  malax  malax-(N/2) ];
        ky = [ malay-(N/2)  malay-(N/2);  malay  malay ];
        ky = flipud(ky);
        kx(:,(N/2+1)) =  kx(:,(N/2+1));
        ky((N/2+1),:) =  ky((N/2+1),:);

        k_abs = sqrt(kx.^2 + ky.^2);
        %%
        alpha = atan2(ky, kx);
        if kx(1,1) == 0 && ky(1,1) == 0,
            alpha(1,1) = 1.57080;
        end;

        %%
        Cx = ((k_abs * young_e12) / (2 * (1 - pois^2))) .* (1 - pois + pois .* (cos(alpha)).^2);
        Cy = ((k_abs * young_e12) / (2 * (1 - pois^2))) .* (1 - pois + pois .* (sin(alpha)).^2);
        D  = ((k_abs * young_e12) / (2 * (1 - pois^2))) .* (pois .* sin(alpha) .* cos(alpha));

        D(:,(N/2+1)) = zeros(N,1);
        D((N/2+1),:) = zeros(1,N);
        %%
        Dx = fft2(ux * 2 * pi / (N * spacing));
        Dy = fft2(uy * 2 * pi / (N * spacing));

        Tx = Cx .* Dx  +  D  .* Dy;
        Ty = D  .* Dx  +  Cy .* Dy;

        tx = real(ifft2(Tx));
        ty = real(ifft2(Ty));

        U1=0.5 * sum(sum([tx ty] .* [ux uy])) * spacing^2 * 1e-6*1e12;

        u_bar=(tx.*ux+ty.*uy)*1e-6;
      U2=0.5*trapz(Y2*1e-6,trapz(X2*1e-6,u_bar))*1e12;

      %interior = find(x2 ~= min(x2(:)) & x2 ~= max(x2(:)) & y2 ~= min(y2(:)) & y2 ~= max(y2(:)));
    u_bar_int=(tx(2:length(tx)-1,2:length(tx)-1).*ux(2:length(tx)-1,2:length(tx)-1)+ty(2:length(tx)-1,2:length(tx)-1).*uy(2:length(tx)-1,2:length(tx)-1))*1e-6;
    Y22=Y2(Y2 ~= Y2(1) & Y2 ~= Y2(length(Y2)));
    X22=X2(X2 ~= X2(1) & X2 ~= X2(length(X2)));

     U3(kk)=0.5*trapz(Y22*1e-6,trapz(X22*1e-6,u_bar_int))*1e12;
     U4=U3.';
     Strainenergy=U3;%Final Strain Energy
     %kk=kk+1;
end
 