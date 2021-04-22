%%%%%% Generates Pauli matrices, Gamma matrices, and spin-j angular momentum operators %%%%%%%%%%%%%%%

%Represent quantum number j by q, m by r
% j is reserved for sqrt(-1), m is reserved for L_x
if ~exist('q','var')
	q = 1;
end
Jz = zeros(2*q + 1, 2*q + 1);
J_minus = zeros(2*q + 1, 2*q +1);
r = q;
for a=1:(2*q)
	Jz(a,a) = r;
	J_minus(a+1,a) = sqrt(q*(q+1) - r*(r-1));
	r = r-1;
endfor
Jz(2*q+1,2*q+1) = -q;
J_plus = transpose(J_minus);

Jx = (J_plus + J_minus)/2;
Jy = (J_plus - J_minus)/(2*j);
antidiag = flipud(eye(2*q + 1));
clear r; %This was only required to generate the matrices

%Pauli matrices
sigma_x = [0, 1; 1, 0];
sigma_y = [0, -j; j, 0];
sigma_z = [1, 0; 0, -1];
sigma_0 = eye(2);
%gamma matrices
gamma_1 = kron(sigma_z, sigma_x);
gamma_2 = kron(sigma_0,sigma_y);
gamma_3 = kron(sigma_0, sigma_z);
gamma_4 = kron(sigma_x, sigma_x);
gamma_5 = kron(sigma_y, sigma_x);