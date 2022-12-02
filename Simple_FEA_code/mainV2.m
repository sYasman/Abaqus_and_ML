clear, clc

% Material Properties
E = 3e7;
nu = 0.30;

% plane stress problem
factor = E / ((1+nu)*(1-2*nu));
factor2 = [1-nu, nu, nu, 0, 0, 0
    nu, 1-nu, nu, 0, 0, 0
    nu, nu, 1-nu, 0, 0, 0
    0, 0, 0, (1-2*nu)/2, 0, 0
    0, 0, 0, 0, (1-2*nu)/2, 0
    0, 0, 0, 0, 0, (1-2*nu)/2];
D_el=factor*factor2;

% read mesh, connectivity and BC data
fid=fopen('n_cords.txt');
nodeCords =cell2mat(textscan(fid,'%f%f%f%f','headerlines',2,'whitespace',', ','collectoutput',1));
fclose(fid);

max_cord = max(nodeCords);

fid=fopen('el_nodes.txt','r');
elNodes =cell2mat(textscan(fid,'%d%d%d%d%d%d%d%d%d','headerlines',2,'whitespace',', ','collectoutput',1));
fclose(fid);

fid=fopen('restrained_nodes.txt');
BC = cell2mat(textscan(fid,'%f%f%f%f%f%f%f','headerlines',2,'whitespace',', ','collectoutput',1));
fclose(fid);

% node and element numbers
n_node = size(nodeCords,1);
n_el = size(elNodes,1);

% update DOF numbering
DOF = zeros(n_node,3);
k = 1;

% number free dof first
for i=1:n_node
    bcflag = 0;
    for j=1:size(BC,1)
        if(i ==  BC(j,1))
            bcflag = 1;
            row_id = j;
        end
    end
    if(bcflag == 0)
        DOF(i,1) = k;
        DOF(i,2) = k+1;
        DOF(i,3) = k+2;
        k = k+3;
    end
    if(bcflag == 1 && BC(row_id,2) == 0)
        DOF(i,1) = k;
        k = k+1;
    end
    if(bcflag == 1 && BC(row_id,3) == 0)
        DOF(i,2) = k;
        k = k+1;
    end
    if(bcflag == 1 && BC(row_id,4) == 0)
        DOF(i,3) = k;
        k = k+1;        
    end
end

% number of not prescribed dof
nfrdof = max(max(DOF));

% numbering prescribed dof
for i=1:n_node
    for j=1:size(DOF,2)
        if(DOF(i,j) == 0)
            DOF(i,j) = k;
            k = k+1;
        end
    end
end

% DOF number
n_dof = max(max(DOF));

% Stiffness matrix
K = zeros(n_dof,n_dof);

% x, y and z coordinates for element
x = zeros(8,1);
y = zeros(8,1);
z = zeros(8,1);

% number of Gauss Points
n_gp = 3;

% element's DOF numbers
C_el = [];

for i=1:n_el % loop over elements  
    elNodes_i = elNodes(i,2:end); % ith elements nodes  
    
    % nodal coordinates for ith elements
    for j=1:size(elNodes_i,2)
        x(j) = nodeCords(elNodes_i(j),2);
        y(j) = nodeCords(elNodes_i(j),3);
        z(j) = nodeCords(elNodes_i(j),4);
    end
    
    % element stiffness matrix
    K_el=getK_el( n_gp, x, y, z, D_el );
    
    % ith elements DOFs
    for k=1:size(elNodes_i,2)
        C_el = [C_el DOF(elNodes_i(k),:)];
    end
    
    % total stiffness matrix
    for j=1:size(C_el,2)
        for k=1:size(C_el,2)
            K(C_el(j),C_el(k)) = K(C_el(j),C_el(k)) + K_el(j,k);
        end
    end
    % flush element DOF array
    C_el = [];
end

% Displacement Vector
U = zeros(n_dof,1);

for i=1:size(BC,1)
    node_w_disp = BC(i,1);
    if(BC(i,2) ~= 0)
        relatedDof = DOF(BC(i,1),1);
        U(relatedDof) = BC(i,5);
    end
    if(BC(i,3) ~= 0)
        relatedDof = DOF(BC(i,1),2);
        U(relatedDof) = BC(i,6);
    end
    if(BC(i,4) ~= 0)
        relatedDof = DOF(BC(i,1),3);
        U(relatedDof) = BC(i,7);
    end    
end

% Prescribed Displacements
Ub = U(nfrdof+1:end);

% split K
% a free DOF
% b DOF with prescribed displacement
Kaa = K((1:nfrdof),(1:nfrdof));
Kab = K((1:nfrdof),(nfrdof+1:end));
Kba = Kab';
Kbb = K((nfrdof+1:end),(nfrdof+1:end));

Ua = Kaa\(-Kab*Ub); % Displacements

Ra = zeros(nfrdof,1);
Rb = Kba*Ua+Kbb*Ub;

U = [Ua; Ub];
node2_m = U(DOF(2,:));
node2_A = [-0.0159103367477655; -0.0417640879750252; -0.000271203345619142];

err = abs(node2_m-node2_A)

