close all
clear all
home

x	   = linspace(0,6,1000);
theta1 = 2;
theta2 = 4;
colors = pgm_colors();

P	   = ( exp(-x/theta1)/theta1.*(x>=0) )' * ( 1/theta2*(x>=0).*(x<theta2) );

%%
figure
hold on; grid
contour(x,x,P,20)
colorbar
plot(1,1,'.','color',colors{5},'MarkerSize',25)
plot(3,3,'.','color',colors{1},'MarkerSize',25)
line([2 2],[0 6],'LineWidth',3, 'color', colors{7})
axis tight
title('\fontsize{16}Curvas de nivel de p(x_1,x_2; \theta_0) y p(x_1,x_2; \theta_1)')
legend('\fontsize{16}Curvas nivel','\fontsize{16}(1,1)^T','\fontsize{16}(3,3)^T','\fontsize{16}(2,*)^T', 'location','northwest', 'orientation', 'horizontal');
xlabel('\fontsize{16}x_1'); ylabel('\fontsize{16}x_2')


%%
[X,Y] = meshgrid(x,x);
P = ( exp(-X/theta1)/theta1.*(X>=0) ).*( 1/theta2*(Y>=0).*(Y<theta2));
figure;
contour(X,Y,P)
