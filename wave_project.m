%syms x y z;
%Inputs
lambda = 3.4e-6; % defined per group (Group 8)
wo = 30e-6; % defined per section (section 3)
k = 2*pi/lambda;
A = 1;
zo = pi*(wo^2)/lambda;
step = sqrt(2)*pi/k;
%x=linspace(-zo*5,zo*5,1000);
%y=linspace(-zo*5,zo*5,1000);
%z=linspace(-zo*1.5,zo*1.5,1000);
x = -wo*5:step:wo*5;
y = -wo*5:step:wo*5;
%Calculations

%k = [kx;ky;kz];
%r = [x;y;z];
%E = @(var) real(A.*exp(-1i.*var));
p = @(x,y) (x.^2+y.^2);
U = @(x,y) A.*exp(-p(x,y)/(wo^2));
%z = 0*wo;
[X,Y] = meshgrid(x,y);
%------------------------------------------
%Part 1
%------------------------------------------
figure(1)
imagesc(x,y,U(X,Y))
caxis([0 1]);
colorbar
colormap jet
%------------------------------------------
%Part 2
%------------------------------------------
figure(2)
subplot(2,1,1);
L = numel(x);
Ts = mean(diff(x));                                            
Fs = 1/step;                                                     
Fn = Fs/2;
xmax=mean(x);
ymax=mean(y);

Ukxky = fftshift(fft2(U(X,Y)));
%[FX,FY] = meshgrid(Fv,Fv);
kx = (-(length(Ukxky)/2):(length(Ukxky)/2)-1).*Fs./length(Ukxky);
kx=kx.*2*pi;
ky = (-(length(Ukxky)/2):(length(Ukxky)/2)-1).*Fs./length(Ukxky);
ky=ky.*2*pi;
imagesc(kx,ky,abs(Ukxky))
colorbar
caxis([0 1]);
subplot(2,1,2);
imagesc(kx,ky,angle(Ukxky))
colorbar
caxis([0 1]);
figure(3)
subplot(2,2,1);
[kxx,kyy] = meshgrid(kx,ky);
% kz = sqrt(k^2.-kxx.^2-kyy.^2);
kz = k - ((kxx.^2+kyy.^2)./(2*k));
U_z = @(z) abs(ifft2(ifftshift(Ukxky.*exp(-1i.*kz.*z))));
imagesc(x,y,(U_z(0.5*zo)).^2)

caxis([0 1]);
colorbar
subplot(2,2,2);
imagesc(x,y,(U_z(1*zo)).^2)
caxis([0 1]);
colorbar
subplot(2,2,3);
imagesc(x,y,(U_z(2*zo)).^2)
caxis([0 1]);
colorbar
colormap jet