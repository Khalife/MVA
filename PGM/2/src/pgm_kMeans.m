function [finalMu, finalLabel] = pgm_kMeans(x, c, opt)
% -------------------------------------------------------------------------
% function [finalMu, finalLabel] = pgm_kMeans(x, c, opt)
% -------------------------------------------------------------------------
% inputs: 
%       - x:   data dxN to be classified. d is the dimension, and N the
%              number of sampels
%       - c:   number of clusters to use
%       - opt: struct with options
% outputs:
%       - finalMu:    matrix of dimensions dxc with the final position of 
%                     the means of each cluster
%       - finalLabel: Nx1 vector of final labels of x
% -------------------------------------------------------------------------

%% Data and initialization

if nargin<3
    opt.plot = 0;
end

[d,N] = size(x);
mu(:,:,1) = 10*randn(d,c); % Randomly initialized means

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

%% Plot
if opt.plot >= 1
    colors = pgm_colors();
    lut    = [1 5 3 10];
    figure
        hold on; grid on;
        for i=1:c
            plot(x(1,l==i),x(2,l==i),'+','color',colors{lut(i)}, 'markersize',5)
            plot(squeeze(mu(1,i,:)),squeeze(mu(2,i,:)),'-^', 'color', 'k', 'MarkerEdgeColor', 'k', 'MarkerFaceColor',colors{lut(i)}, 'markersize', 8, 'linewidth', 1)
        end
    axis tight
end
finalMu = mu(:,:,end);
finalLabel = l;