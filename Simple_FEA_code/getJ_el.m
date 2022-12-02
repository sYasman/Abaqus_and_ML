function [ Jacobian, B_strain_disp ] = getJ_el( x, y, z, xi, eta, tau )
% Function returns Jacobina for linear 4 point 2D elements
% x, y, z element node coordinates ! Must be column vector
% xi, eta, tau gauss points

%N_grad is grad(xi, eta, tau).N
N_grad = zeros(3,8);
N_grad(1,1) = -(1/8)*(1-eta)*(1+tau); % -(1/8) (1-eta) (1+tau)
N_grad(1,2) = (1/8)*(1-eta)*(1+tau); % 1/8 (1-eta) (1+tau)
N_grad(1,3) = (1/8)*(1-eta)*(1-tau); % 1/8 (1-eta) (1-tau)
N_grad(1,4) = -(1/8)*(1-eta)*(1-tau); % -(1/8) (1-eta) (1-tau
N_grad(1,5) = -(1/8)*(1+eta)*(1+tau); % -(1/8) (1+eta) (1+tau)
N_grad(1,6) = (1/8)*(1+eta)*(1+tau); % 1/8 (1+eta) (1+tau)
N_grad(1,7) = (1/8)*(1+eta)*(1-tau); % 1/8 (1+eta) (1-tau)
N_grad(1,8) = -(1/8)*(1+eta)*(1-tau); % -(1/8) (1+eta) (1-tau)

N_grad(2,1) = -(1/8)*(1-xi)*(1+tau); % -(1/8) (1-xi) (1+tau)
N_grad(2,2) = -(1/8)*(1+xi)*(1+tau); % -(1/8) (1+xi) (1+tau)
N_grad(2,3) = -(1/8)*(1+xi)*(1-tau); % -(1/8) (1+xi) (1-tau)
N_grad(2,4) = -(1/8)*(1-xi)*(1-tau); % -(1/8) (1-xi) (1-tau)
N_grad(2,5) = (1/8)*(1-xi)*(1+tau); % 1/8 (1-xi) (1+tau)
N_grad(2,6) = (1/8)*(1+xi)*(1+tau); % 1/8 (1+xi) (1+tau)
N_grad(2,7) = (1/8)*(1+xi)*(1-tau); % 1/8 (1+xi) (1-tau
N_grad(2,8) = (1/8)*(1-xi)*(1-tau); % 1/8 (1-xi) (1-tau)

N_grad(3,1) = (1/8)*(1-eta)*(1-xi); % 1/8 (1-eta) (1-xi)
N_grad(3,2) = (1/8)*(1-eta)*(1+xi); % 1/8 (1-eta) (1+xi)
N_grad(3,3) = -(1/8)*(1-eta)*(1+xi); % -(1/8) (1-eta) (1+xi)
N_grad(3,4) = -(1/8)*(1-eta)*(1-xi); % -(1/8) (1-eta) (1-xi)
N_grad(3,5) = (1/8)*(1+eta)*(1-xi); % 1/8 (1+eta) (1-xi)
N_grad(3,6) = (1/8)*(1+eta)*(1+xi); % 1/8 (1+eta) (1+xi)
N_grad(3,7) = -(1/8)*(1+eta)*(1+xi); % -(1/8) (1+eta) (1+xi)
N_grad(3,8) = -(1/8)*(1+eta)*(1-xi); % -(1/8) (1+eta) (1-xi)

% Jacobian matrix
Jacobian = N_grad * [x,y,z];

% B matrix : [Ni,x Ni,y, Ni,z]' = J^-1 . Grad(N, (xi, eta))
B_part = Jacobian\N_grad;
B_strain_disp = zeros(6,24);

N_x = B_part(1,:);
N_y = B_part(2,:);
N_z = B_part(3,:);
% k=1;
% for ith element, i increases by 3 up to 22
% B : N_x(1,i),    0,          0
%     0,        N_y(1,i),      0
%     0,           0,      N_z(1,i)
%     N_y(1,i), N_x(1,i),      0
%     N_z(1,i),    0,      N_x(1,i)
%     0,        N_z(1,i),  N_y(1,i)

% for i=1:3:22
%     B(1:6,i:i+2) = [N_x(k), 0 , 0
%         0, N_y(k), 0
%         0, 0, N_z(k)
%         N_y(k), N_x(k), 0
%         N_z(k), 0, N_x(k)
%         0, N_z(k), N_y(k)];
%     k=k+1;
% end

B_strain_disp(1,1)=N_x(1);
B_strain_disp(1,4)=N_x(2);
B_strain_disp(1,7)=N_x(3);
B_strain_disp(1,10)=N_x(4);
B_strain_disp(1,13)=N_x(5);
B_strain_disp(1,16)=N_x(6);
B_strain_disp(1,19)=N_x(7);
B_strain_disp(1,22)=N_x(8);

B_strain_disp(2,2)=N_y(1);
B_strain_disp(2,5)=N_y(2);
B_strain_disp(2,8)=N_y(3);
B_strain_disp(2,11)=N_y(4);
B_strain_disp(2,14)=N_y(5);
B_strain_disp(2,17)=N_y(6);
B_strain_disp(2,20)=N_y(7);
B_strain_disp(2,23)=N_y(8);

B_strain_disp(3,3)=N_z(1);
B_strain_disp(3,6)=N_z(2);
B_strain_disp(3,9)=N_z(3);
B_strain_disp(3,12)=N_z(4);
B_strain_disp(3,15)=N_z(5);
B_strain_disp(3,18)=N_z(6);
B_strain_disp(3,21)=N_z(7);
B_strain_disp(3,24)=N_z(8);

B_strain_disp(4,1)=N_y(1);
B_strain_disp(4,4)=N_y(2);
B_strain_disp(4,7)=N_y(3);
B_strain_disp(4,10)=N_y(4);
B_strain_disp(4,13)=N_y(5);
B_strain_disp(4,16)=N_y(6);
B_strain_disp(4,19)=N_y(7);
B_strain_disp(4,22)=N_y(8);

B_strain_disp(4,2)=N_x(1);
B_strain_disp(4,5)=N_x(2);
B_strain_disp(4,8)=N_x(3);
B_strain_disp(4,11)=N_x(4);
B_strain_disp(4,14)=N_x(5);
B_strain_disp(4,17)=N_x(6);
B_strain_disp(4,20)=N_x(7);
B_strain_disp(4,23)=N_x(8);

B_strain_disp(5,1)=N_z(1);
B_strain_disp(5,4)=N_z(2);
B_strain_disp(5,7)=N_z(3);
B_strain_disp(5,10)=N_z(4);
B_strain_disp(5,13)=N_z(5);
B_strain_disp(5,16)=N_z(6);
B_strain_disp(5,19)=N_z(7);
B_strain_disp(5,22)=N_z(8);

B_strain_disp(5,3)=N_x(1);
B_strain_disp(5,6)=N_x(2);
B_strain_disp(5,9)=N_x(3);
B_strain_disp(5,12)=N_x(4);
B_strain_disp(5,15)=N_x(5);
B_strain_disp(5,18)=N_x(6);
B_strain_disp(5,21)=N_x(7);
B_strain_disp(5,24)=N_x(8);

B_strain_disp(6,2)=N_z(1);
B_strain_disp(6,5)=N_z(2);
B_strain_disp(6,8)=N_z(3);
B_strain_disp(6,11)=N_z(4);
B_strain_disp(6,14)=N_z(5);
B_strain_disp(6,17)=N_z(6);
B_strain_disp(6,20)=N_z(7);
B_strain_disp(6,23)=N_z(8);

B_strain_disp(6,3)=N_y(1);
B_strain_disp(6,6)=N_y(2);
B_strain_disp(6,9)=N_y(3);
B_strain_disp(6,12)=N_y(4);
B_strain_disp(6,15)=N_y(5);
B_strain_disp(6,18)=N_y(6);
B_strain_disp(6,21)=N_y(7);
B_strain_disp(6,24)=N_y(8);

end

