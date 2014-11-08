close all
clear all
home

%% Data and initialization
opt.plot = 1;
opt.log  = 1;

x        = load('../data/EMGaussian.data'); x = x';
[d,N]    = size(x);
c        = 4; % Number of clusters
if opt.log>=1
    fprintf('Running K-Means...\n');
end
mu       = pgm_kMeans(x, c, opt)';
if opt.log>=1
    fprintf('\tDone!\n');
end
t        = 1;
epsilon  = 1e-03;
Pi       = (1/c)*ones(c,1);
SIGMA    = zeros(d,d,c);

for j=1:c
    SIGMA(:,:,j)=eye(d);
end

%% Expectation maximization
if opt.log>=1
    fprintf('Running Expectation Maximization algorithm...\n');
end
l        = -inf;
lDiff    = inf;
counter  = 1;
while (lDiff > epsilon)
    % Given theta, with bayes formula compute p(z|x)=q* (partial maximisation)
    p_z_x = pgm_compute_p_z_x(x,Pi,mu,SIGMA);   
    l_new = pgm_computeLikelyhood(p_z_x,x,Pi,mu,SIGMA)
    lDiff = abs(l-l_new);
    % Closed formula for updating the theta variables given p(z|x)
    Pi    = sum(p_z_x)/N;  
    % Mean actualization
    for j=1:c
        mu(j,:) = sum(repmat(p_z_x(:,j),1,2).*x')./sum(p_z_x(:,j));
    end
    % Sigma actualization
    for j=1:c
        SIGMA(:,:,j)=zeros(d,d);
        for i=1:N
            SIGMA(:,:,j) = SIGMA(:,:,j) + (p_z_x(i,j))*(x(:,i)'-mu(j,:))'*(x(:,i)'-mu(j,:));
        end
        SIGMA(:,:,j) = SIGMA(:,:,j)./sum(p_z_x(:,j));
    end
    if opt.log>=2
        fprintf('Iteration %d:\tliklyhood %0.4f\n', counter, l_new);
    end
    l       = l_new;
    counter = counter+1;
end
if opt.log>=1
    fprintf('\tDone!\n');
end

%% Plot

colors = pgm_colors();
lut    = [1 5 3 10];
    
if opt.log>=1
    fprintf('Plotting...\n');
end

if opt.plot>=2
    pdf   = @(x,mu,SIGMA) exp(-0.5*(x-mu)'*(SIGMA\(x-mu)))/sqrt(det(SIGMA)*(2*pi)^d);
    xx = linspace(min(x(1,:)),max(x(1,:)),100);
    yy = linspace(min(x(2,:)),max(x(2,:)),100);
    [X,Y] = meshgrid(xx,yy);
    surface1 = zeros(size(X));
    surface2 = zeros(size(X));
    surface3 = zeros(size(X));
    surface4 = zeros(size(X));
    for i=1:size(X,1)
        for j=1:size(X,2)
            surface1(i,j) = pdf([X(i,j);Y(i,j)],mu(1,:)',SIGMA(:,:,1));
            surface2(i,j) = pdf([X(i,j);Y(i,j)],mu(2,:)',SIGMA(:,:,2));
            surface3(i,j) = pdf([X(i,j);Y(i,j)],mu(3,:)',SIGMA(:,:,3));
            surface4(i,j) = pdf([X(i,j);Y(i,j)],mu(4,:)',SIGMA(:,:,4));
        end
    end
    figure;
    plot(x(1,:),x(2,:),'+')
    hold on
    contour(X,Y,surface1)
    contour(X,Y,surface2)
    contour(X,Y,surface3)
    contour(X,Y,surface4)
    grid on
    axis([xx(1) xx(end) yy(1) yy(end)]);
    title('LEVEL SETS')
end
if opt.plot>=1
    figure;
    plot(x(1,:),x(2,:),'+','markersize',5)
    hold on
    grid on
for i=1:c
    [xEllip, yEllip] = pgm_computeEllipse(mu(i,:),SIGMA(:,:,i),0.9);
    plot(xEllip, yEllip,'color',colors{lut(i)},'linewidth',2)
    axis tight
    title('\fontsize{14}90% of the mass of the Gaussian')
end
end

if opt.log>=1
    fprintf('\tDone!\n');
end
