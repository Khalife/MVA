function [xEllip, yEllip] = pgm_computeEllipse(mu, SIGMA, prob)

sigma1 = sqrt(SIGMA(1,1));
sigma2 = sqrt(SIGMA(2,2));
rho    = SIGMA(1,2)/sigma1/sigma2;
r      = sqrt(-2*log(1-prob));
theta  = linspace(0,2*pi,1000);

xEllip = r*sigma1*cos(theta)+mu(1);
yEllip = r*sigma2*(rho*cos(theta)+sqrt(1-rho^2)*sin(theta))+mu(2);