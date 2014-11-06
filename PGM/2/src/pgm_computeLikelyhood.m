function l=pgm_computeLikelyhood(p_z_x,x,Pi,mu,SIGMA)

l=0;
for i=1:size(x,2)
    for j=1:size(p_z_x,2)
        l=l+p_z_x(i,j)*(log(Pi(j))-log(2*pi)-0.5*log(det(SIGMA(:,:,j)))-0.5*(x(:,i)'-mu(j,:))*(SIGMA(:,:,j)\(x(:,i)'-mu(j,:))'));
    end
end

end