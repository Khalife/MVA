close all
clear all
home

%% Data and initialization

x     = load('EMGaussian.data'); x = x';
[d,N] = size(x);

c         = 4;          % Number of clusters
mu(:,:) = randn(d,c); % Randomly initialized means

t = 1;
epsilon = 1e-07;
n=size(x,2);
%Initialisation of parameters
Pi=(1/4)*ones(c,1);
mu=zeros(c,2);
SIGMA=zeros(d,d,c);
for j=1:c
SIGMA(:,:,j)=eye(d,d);
end

%p_z_x=compute_p_z_x(x,Pi,mu,SIGMA);
%l=computeLikelyhood(p_z_x,x,Pi,mu,SIGMA);
compteur=1;
while (1)
%given theta, with bayes formula compute p(z|x)=q* (corresponds to a partial maximisation) 
p_z_x=compute_p_z_x(x,Pi,mu,SIGMA);

l_new=computeLikelyhood(p_z_x,x,Pi,mu,SIGMA)
 
if ( compteur > 1 && abs(l-l_new) < epsilon )
    
    break;
end

%closed formula for updating the theta variables given p(z|x)
for j=1:c
Pi(j)=(1/c)*(sum(p_z_x(j,:)));
end

for j=1:c
mu(j,1)=0;
mu(j,2)=0;
    for i=1:N
        mu(j,:) = mu(j,:) + [p_z_x(i,j)*x(1,i) p_z_x(i,j)*x(2,i)]; 
    end
mu(j,:) = mu(j,:)./sum(p_z_x(:,j));
end

for j=1:c
SIGMA(:,:,j)=zeros(d,d);
    
    for i=1:n
        SIGMA(:,:,j) = SIGMA(:,:,j) + (p_z_x(i,j))*(x(:,i)'-mu(j,:))'*(x(:,i)'-mu(j,:));
    end
    
SIGMA(:,:,j)= SIGMA(:,:,j)./sum(p_z_x(:,j));
end

l=l_new;
compteur=compteur+1;
end