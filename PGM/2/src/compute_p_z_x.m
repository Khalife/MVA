function p_z_x=compute_p_z_x(x,Pi,mu,sigma)

%sigma is a vector of covariance matrices
%mu a vector of means
n=size(x,2);
c=size(mu,1);
p_z_x=zeros(n,c);
for i=1:n
    for j=1:c
        p_z_x(i,j)=Pi(j)*Normale(x(:,i),mu(j,:),sigma(:,:,j));
    end
    p_z_x(i,:)=p_z_x(i,:)./sum(p_z_x(i,:));
end

end