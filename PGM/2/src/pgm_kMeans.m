close all
clear all
home

%% Data and initialization

x     = load('../data/EMGaussian.data'); x = x';
[d,N] = size(x);

c         = 4;          % Number of clusters
mu(:,:,1) = randn(d,c); % Randomly initialized means

t = 1;
J = 1;
l = zeros(N,1);
d = zeros(c,N);

%% loop k-means
while J > 0 && t < 1000
	J = 0;
	for i=1:c
		aux    = x-repmat(mu(:,i,end),1,N);
		d(i,:) = aux(1,:).^2+aux(2,:).^2;
	end
	[~,l] = min(d);
	
	for i=1:c
		mu(:,i,t+1) = mean(x(:,l==i),2);
		J = J+(mu(:,i,t+1)-mu(:,i,t)).^2;
	end
	J = norm(J);
    t = t+1;
end

%% Plots
colors = pgm_colors();
lut    = [1 5 3 10];
figure
	hold on; grid on;
	for i=1:c
		plot(x(1,l==i),x(2,l==i),'.','color',colors{lut(i)}, 'markersize',15)
		plot(squeeze(mu(1,i,:)),squeeze(mu(2,i,:)),'-^', 'color', 'k', 'MarkerEdgeColor', 'k', 'MarkerFaceColor',colors{lut(i)}, 'markersize', 10, 'linewidth', 2)
	end
axis tight