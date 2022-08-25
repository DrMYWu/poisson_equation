% Solving 2D Poisson equation By Gauss iterative method.
%All units are in SI.
%Writed by Wu Ming Yang.
% Date:2017/11/13

clear;clc;close all;
qe=1.6e-19;%C
ep0=8.89e-12;
%%
center_rho_f_m3=1e12; %charge density,m-3
%center_ni_m3=0;
delta=1e-4; %precision

%Calculation area
xs=0;xd=4;
ys=0;yd=4;

%Number of grids
nx=100;
ny=100;

x(nx)=0;
y(ny)=0;
dx=(xd-xs)/(nx-1);
dy=(yd-ys)/(ny-1);
for ix=1:nx
    x(ix)=xs+(ix-1)*dx;
end
for iy=1:ny
    y(iy)=ys+(iy-1)*dy;
end

x1=1;y1=1;v1=0;%V   (x1,y1)~(x3,y1) votage is set to v1
x3=3;y3=3;v3=100;%V (x1,y3)~(x3,y3) votage is set to v3
x2=2;y2=2;%position of charge density

nx1=floor((x1-xs)/dx);
nx2=floor((x2-xs)/dx);
nx3=floor((x3-xs)/dx);
ny1=floor((y1-ys)/dy);
ny2=floor((y2-ys)/dy);
ny3=floor((y3-ys)/dy);

v=zeros(nx,ny);
vo=zeros(nx,ny);
rho_f=zeros(nx,ny);

dx2=1/dx^2;
dy2=1/dy^2;
dxy=2*(1/dx^2+1/dy^2);

v(nx1,ny1:ny3)=v1;
v(nx3,ny1:ny3)=v3;
rho_f(nx2,ny2)=qe*center_rho_f_m3/ep0; %charge density

k=0;
while 1
    vtmp=v;
    for i=2:nx-1
        for j=2:ny-1
            tp1=dx2*v(i+1,j)+dx2*v(i-1,j)+dy2*v(i,j+1)+dy2*v(i,j-1)+rho_f(i,j);
            v(i,j)=tp1/dxy;
        end
    end

    %boundary condition
    v(1,1:ny)=v(2,1:ny);
    v(nx,1:ny)=v(nx-1,1:ny);

    v(1:nx,1)=v(1:nx,2);
    v(1:nx,ny)=v(1:nx,ny-1);

    %Voltage setting
    v(nx1,ny1:ny3)=v1;
    v(nx3,ny1:ny3)=v3;

    % condition fo Jumping out of the loop 
    if max(abs(v-vtmp))<delta
        break;
    end
    k=k+1;
end

Number_of_iterations=k
%%
i_switch_print=1;%1-sve the figure; 0-not save
resolution='-r450';
ls=25;lw=3;

ps = get(0, 'ScreenSize');
ps1(1) = floor(ps(3)*0.02);
ps1(2) = floor(ps(4)*0.5);
ps1(3) = floor(ps(3)*0.3);
ps1(4) = floor(ps(4)*0.4);

ps2(1) = floor(ps(3)*0.02);
ps2(2) = floor(ps(4)*0.05);
ps2(3) = floor(ps(3)*0.3);
ps2(4) = floor(ps(4)*0.4);

%%
figure
mesh(x,y,v)
xlabel('x (m)');
ylabel('y (m)');
zlabel('V (V)');
% title('The potential of pallel capacitor');
set(gca, 'fontsize', ls);
set(gcf,'position',ps1)
fname='V_xy_3D';
if i_switch_print>0.1
    print([fname,'.png'],'-dpng',resolution);
end


figure
contour(x,y,v,20)
xlabel('x (m)');
ylabel('y (m)');
zlabel('V (m)');
% title('The potential of pallel capacitor');
set(gca, 'fontsize', ls);
set(gcf,'position',ps2)
fname='V_xy_contour';
if i_switch_print>0.1
    print([fname,'.png'],'-dpng',resolution);
end



