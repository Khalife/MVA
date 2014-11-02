function n=Normale(x,mu,sigma)


d=size(mu,2);
n=(1/((det(sigma))^(0.5)*(2*pi)^(d/2)))*exp(-(1/2)*(x-mu')'*(sigma\(x-mu')));

end