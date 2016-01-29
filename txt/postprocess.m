function [] = postprocess(proc)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
X = load(strcat('process_xcoord_',num2str(proc),'.txt'));
Y = load(strcat('process_ycoord_',num2str(proc),'.txt'));
P = load(strcat('pressure_',num2str(proc),'.txt'));
figure(1)
hold on;
mesh(X,Y,P);
xlim([-3.5 3.5])
ylim([0 2])
colorbar
%title(strcat('Pressure 2D plot for M = ',sprintf('%g',mInfi)),'FontSize', 14);
xlabel('X-Coord');
ylabel('Y-Coord');
figure(2)
hold on
h = pcolor(X,Y,P);
set(h,'EdgeColor','none','LineStyle','none','FaceLighting','phong');
shading interp
xlim([-3.5 3.5])
ylim([0 2])
colorbar
%title(strcat('Pressure 2D plot for M = ',sprintf('%g',mInfi)),'FontSize', 14);
xlabel('X-Coord');
ylabel('Y-Coord');

end

