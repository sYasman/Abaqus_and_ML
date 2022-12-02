function [ K_el ] = getK_el( n_gp, x, y, z, D_el )
% function returns element stiffness matrix by using Gauss Quadrature
% getK_el( n_gp, x, y, z, D_el );
% n_gp : Gauss Points for element
% x : x coordinates of element
% y : y coordinates of element
% D_el : stiffness matrix

K_el = zeros(24,24);
[w, xi_vals] = getGaussData( n_gp );

for i=1:n_gp
    for j=1:n_gp
        for k=1:n_gp
            [J_el, B_el] = getJ_el( x, y, z, xi_vals(i), xi_vals(j), xi_vals(k) );
            K_el = K_el + w(i)*w(j)*w(k)*B_el'*D_el*B_el*det(J_el);
        end
    end
end




end


