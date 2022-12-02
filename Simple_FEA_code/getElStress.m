function [ stress_el ] = getElStress( no_el, elNodes, D_el, DOF, U, n_gp )
% function returns stress at Gauss Point
% no_el element number
% elNodes nodes of elements (el_noi node_1, node_2, node_3, node_4)
% D_el Element stiffness matrix
% DOF DOF matrix (Dof_1 Dof_2)
% U displacement vector
% n_gp Number of Gauss Integration points
no_el = 1;
u_el = [];
el_dof = [];
nodes = elNodes(no_el,2:end);

for k=1:size(nodes,2)
    el_dof = [el_dof DOF(k,:)];
    % ith element x and y displacements
    u_el = [u_el U(DOF(nodes(k),1)) U(DOF(nodes(k),2))];
end

if n_gp == 2
    xi_vals = [1/sqrt(3), -1/sqrt(3)];
    eta_vals = [1/sqrt(3), -1/sqrt(3)];
end
stress_el = zeros(3, n_gp*n_gp);
p1=1;
% strains and stresses at Gauss Points
for ii=1:n_gp
    for jj=1:n_gp
        B_el = getB_el( xi_vals(ii),eta_vals(jj));
        e_el = B_el * u_el';
        stress_el(:,p1) = D_el * e_el;
        p1=p1+1;
    end
end

end

