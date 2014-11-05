function p_z_x=pgm_compute_p_z_x(x,Pi,mu,SIGMA)
% -------------------------------------------------------------------------
% function p_z_x=compute_p_z_x(x,Pi,mu,SIGMA)
% -------------------------------------------------------------------------
% inputs:
%   - x:  matrix dxN of samples
%   - Pi: vector cx1
%   - mu: matrix cxd of means of clusters
% output:
%   - p_z_x
% -------------------------------------------------------------------------

[d,N] = size(x);
c     = size(mu,1);
pdf   = @(x,mu,SIGMA) exp(-0.5*(x-mu)'*(SIGMA\(x-mu)))/sqrt(det(SIGMA)*(2*pi)^d);
p_z_x = zeros(N,c);

for j=1:c
    for i=1:N
        p_z_x(i,j) = Pi(j)*pdf(x(:,i),mu(j,:)',SIGMA(:,:,j));
    end
end
p_z_x = p_z_x ./ repmat(sum(p_z_x,2),1,c);