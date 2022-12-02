function [ B_el ] = getB_el( xi,eta )
% Return B matrix for given xi and eta values
% use it for gauss points of xi_i, and eta_i
B_el=zeros(3,8);
B_el(1,1) = (1-0.50*eta-0.50*xi)/(eta-3.0);
B_el(1,3) = (0.50*xi+0.50)/(eta-3.0);
B_el(1,5) = (xi+1)/(3-eta);
B_el(1,7) = (0.50-0.50*eta-xi)/(3-eta);

B_el(2,2) = 2 - 4 /(3-eta);
B_el(2,4) = -2 + 4 /(3-eta);
B_el(2,6) = 2 - 8/(3-eta);
B_el(2,8) = -2 + 8 /(3-eta);

B_el(3,1) = 2 - 4/(3-eta);
B_el(3,2) = (1-0.50*eta-0.50*xi)/(eta-3);
B_el(3,3) = -2 + 4/(3-eta);
B_el(3,4) = (0.50+0.50*xi)/(eta-3);
B_el(3,5) = 2 -8/(3-eta);
B_el(3,6) = (1+xi)/(3-eta);
B_el(3,7) = -2 + 8/(3-eta);
B_el(3,8) = (0.50-0.50*eta-xi)/(3-eta);

end

