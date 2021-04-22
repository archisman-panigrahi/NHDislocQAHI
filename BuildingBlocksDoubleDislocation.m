%%%% Generates sine and cosine matrices for a lattice with two dislocation centers %%%%%

%Show lattice size

Lx,Ly,l1,l2,m

n_atoms = Lx * Ly + l2;

% l denotes the row in which dislocation is placed.
% l belongs to the set (1,2,... Ly -1)
% l = 0 means dislocation leaves the system through the bottom layer
% l = Ly means the dislocation leaves the system through the top layer

SXnp = zeros(n_atoms,n_atoms);
SYnp = zeros(n_atoms,n_atoms);
CXnp = zeros(n_atoms,n_atoms);
CYnp = zeros(n_atoms,n_atoms);
% SXYnp = zeros(n_atoms,n_atoms);

M2D = eye(n_atoms);

%%%%%%%%% kx matrices %%%%%%%%%%%
%kx matrices, lower block
for a = 1:(Lx - 1)
	for b = 0:(l1-1)
		SXnp(a + b*Lx, a + b*Lx + 1) = j/2;
		SXnp(a + b*Lx + 1, a + b*Lx) = -j/2;
		CXnp(a + b*Lx, a + b*Lx + 1) = 1/2;
		CXnp(a + b*Lx + 1, a + b*Lx) = 1/2;
	endfor
endfor

%kx matrices, middle block
for a = 1:Lx
	for b = 0:(l2 - 1)
		SXnp(a + l1*Lx + b * (Lx + 1), a + l1*Lx + b * (Lx + 1) + 1) = j/2;
		SXnp(a + l1*Lx + b * (Lx + 1) + 1, a + l1*Lx + b * (Lx + 1)) = -j/2;
		CXnp(a + l1*Lx + b * (Lx + 1), a + l1*Lx + b * (Lx + 1) + 1) = 1/2;
		CXnp(a + l1*Lx + b * (Lx + 1) + 1, a + l1*Lx + b * (Lx + 1)) = 1/2;
	endfor
endfor

%kx matrices, upper block
for a = 1:(Lx-1)
	for b = 0:(Ly - (l1 + l2) - 1)
		SXnp(a + (l1 + l2)*Lx + l2 + b * Lx, a + (l1 + l2)*Lx + l2 + b * Lx + 1) = j/2;
		SXnp(a + (l1 + l2)*Lx + l2 + b * Lx + 1, a + (l1 + l2)*Lx + l2 + b * Lx) = -j/2;
		CXnp(a + (l1 + l2)*Lx + l2 + b * Lx, a + (l1 + l2)*Lx + l2 + b * Lx + 1) = 1/2;
		CXnp(a + (l1 + l2)*Lx + l2 + b * Lx + 1, a + (l1 + l2)*Lx + l2 + b * Lx) = 1/2;
	endfor
endfor

%%% Periodic Boundary Condition %%%
SXp = SXnp;
CXp = CXnp;
%lower block
for b = 0:(l1 - 1)
	SXp((b+1)*Lx, b*Lx + 1) = j/2;
	SXp(b*Lx + 1, (b+1)*Lx) = -j/2;
	CXp((b+1)*Lx, b*Lx + 1) = 1/2;
	CXp(b*Lx + 1, (b+1)*Lx) = 1/2;
endfor

%middle block
for b = 0:(l2 - 1)
	SXp(l1*Lx + (b+1)*(Lx+1), l1*Lx + b*(Lx+1) + 1) = j/2;
	SXp(l1*Lx + b*(Lx+1) + 1, l1*Lx + (b+1)*(Lx+1)) = -j/2;
	CXp(l1*Lx + (b+1)*(Lx+1), l1*Lx + b*(Lx+1) + 1) = 1/2;
	CXp(l1*Lx + b*(Lx+1) + 1, l1*Lx + (b+1)*(Lx+1)) = 1/2;
endfor

%upper block
for b = 0:(Ly - (l1 + l2) -1)
	SXp((b+1)*Lx + Lx*(l1+l2) + l2, b*Lx + 1 + Lx*(l1+l2) + l2) = j/2;
	SXp(b*Lx + 1 + Lx*(l1+l2) + l2, (b+1)*Lx + Lx*(l1+l2) + l2) = -j/2;
	CXp((b+1)*Lx + Lx*(l1+l2) + l2, b*Lx + 1 + Lx*(l1+l2) + l2) = 1/2;
	CXp(b*Lx + 1 + Lx*(l1+l2) + l2, (b+1)*Lx + Lx*(l1+l2) + l2) = 1/2;
endfor
%%%%%%%%%% ky matrices %%%%%%%%%%%
%ky matrices, lower block (without connectors)
for a = 1:Lx
	for b = 0:(l1-2)
		SYnp(a + b*Lx, a + (b+1)*Lx) = j/2;
		SYnp(a + (b+1)*Lx, a + b*Lx) = -j/2;
		CYnp(a + b*Lx, a + (b+1)*Lx) = 1/2;
		CYnp(a + (b+1)*Lx, a + b*Lx) = 1/2;
	endfor
endfor

%middle block, without connectors
for a = 1:(Lx+1)
	for b = 0:(l2-2)
		SYnp(a + b*(Lx+1) + l1*Lx, a + (b+1)*(Lx+1) + l1*Lx) = j/2;
		SYnp(a + (b+1)*(Lx+1) + l1*Lx, a + b*(Lx+1) + l1*Lx) = -j/2;
		CYnp(a + b*(Lx+1) + l1*Lx, a + (b+1)*(Lx+1) + l1*Lx) = 1/2;
		CYnp(a + (b+1)*(Lx+1) + l1*Lx, a + b*(Lx+1) + l1*Lx) = 1/2;
	endfor
endfor

%ky matrices, upper block (without connectors)
for a = 1:Lx
	for b = 0:(Ly - (l1 + l2) - 2)
		SYnp(a + b*Lx + (l1+l2)*Lx + l2, a + (b+1)*Lx + (l1+l2)*Lx + l2) = j/2;
		SYnp(a + (b+1)*Lx + (l1+l2)*Lx + l2, a + b*Lx + (l1+l2)*Lx + l2) = -j/2;
		CYnp(a + b*Lx + (l1+l2)*Lx + l2, a + (b+1)*Lx + (l1+l2)*Lx + l2) = 1/2;
		CYnp(a + (b+1)*Lx + (l1+l2)*Lx + l2, a + b*Lx + (l1+l2)*Lx + l2) = 1/2;
	endfor
endfor

%ky matrices, lower connectors
for a = 1:m
	SYnp(a + (l1-1)*Lx, a + l1*Lx) = j/2;
	SYnp(a + l1*Lx, a + (l1-1)*Lx) = -j/2;
	CYnp(a + (l1-1)*Lx, a + l1*Lx) = 1/2;
	CYnp(a + l1*Lx, a + (l1-1)*Lx) = 1/2;
endfor

for a = 1:(Lx - m)
	SYnp(a + m + (l1-1)*Lx, a + m + 1 + l1*Lx) = j/2;
	SYnp(a + m + 1 + l1*Lx, a + m + (l1-1)*Lx) = -j/2;
	CYnp(a + m + (l1-1)*Lx, a + m + 1 + l1*Lx) = 1/2;
	CYnp(a + m + 1 + l1*Lx, a + m + (l1-1)*Lx) = 1/2;
endfor
%ky matrices, upper connectors
for a = 1:m
	SYnp(a + (l1 + l2)*Lx + l2 - (Lx + 1), a + (l1 + l2)*Lx + l2) = j/2;
	SYnp(a + (l1 + l2)*Lx + l2, a + (l1 + l2)*Lx + l2 - (Lx + 1)) = -j/2;
	CYnp(a + (l1 + l2)*Lx + l2 - (Lx + 1), a + (l1 + l2)*Lx + l2) = 1/2;
	CYnp(a + (l1 + l2)*Lx + l2, a + (l1 + l2)*Lx + l2 - (Lx + 1)) = 1/2;
endfor

for a = 1:(Lx - m)
	SYnp(a + (m+1) + (l1 + l2)*Lx + l2 - (Lx + 1), a + m + (l1 + l2)*Lx + l2) = j/2;
	SYnp(a + m + (l1 + l2)*Lx + l2, a + (m+1) + (l1 + l2)*Lx + l2 - (Lx + 1)) = -j/2;
	CYnp(a + (m+1) + (l1 + l2)*Lx + l2 - (Lx + 1), a + m + (l1 + l2)*Lx + l2) = 1/2;
	CYnp(a + m + (l1 + l2)*Lx + l2, a + (m+1) + (l1 + l2)*Lx + l2 - (Lx + 1)) = 1/2;
endfor

%%%%Periodic Boundary condition%%%%%
SYp = SYnp;
CYp = CYnp;
for a = 1:Lx
	SYp(Lx*(Ly-1) + l2 + a, a) = j/2;
	SYp(a, Lx*(Ly-1) + l2 + a) = -j/2;
	CYp(Lx*(Ly-1) + l2 + a, a) = 1/2;
	CYp(a, Lx*(Ly-1) + l2 + a) = 1/2;
endfor

% %sin(kx)*sin(ky) matrices
% %lower part
% for a = 1:(Lx - 1)
% 	for b = 0:(l1-2)
% 		SXYnp(a + b*Lx, (a+1) + (b+1)*Lx) = -1/4;
% 		SXYnp((a+1) + (b+1)*Lx, a + b*Lx) = -1/4;

% 		SXYnp(a + 1 + b*Lx, a + (b+1)*Lx) = 1/4;
% 		SXYnp(a + (b+1)*Lx, a + 1 + b*Lx) = 1/4;
% 	endfor
% endfor

% %Middle part
% for a = 1:Lx
% 	for b = 0:(l2-2)
% 		SXYnp(a + b*(Lx+1) + l1*Lx, (a+1) + (b+1)*(Lx+1) + l1*Lx) = -1/4;
% 		SXYnp((a+1) + (b+1)*(Lx+1) + l1*Lx, a + b*(Lx+1) + l1*Lx) = -1/4;

% 		SXYnp(a + 1 + b*(Lx+1) + l1*Lx, a + (b+1)*(Lx+1) + l1*Lx) = 1/4;
% 		SXYnp(a + (b+1)*(Lx+1) + l1*Lx, a + 1 + b*(Lx+1) + l1*Lx) = 1/4;
% 	endfor
% endfor

% %Upper part
% for a = 1:(Lx-1)
% 	for b = 0:(Ly - (l1 + l2) -2)
% 		SXYnp(a + b*Lx + (l1 + l2)*Lx + l2, (a+1) + (b+1) * Lx + (l1 + l2)*Lx + l2) = -1/4;
% 		SXYnp((a+1) + (b+1)*Lx + (l1 + l2)*Lx + l2, a + b*Lx + (l1 + l2)*Lx + l2) = -1/4;

% 		SXYnp(a + 1 + b*Lx + (l1+l2)*Lx + l2, a + (b+1) * Lx + (l1+l2)*Lx + l2) = 1/4;
% 		SXYnp(a + (b + 1) * Lx + (l1+l2)*Lx + l2, a + 1 + b*Lx + (l1+l2)*Lx + l2) = 1/4;
% 	endfor
% endfor

% %Connectors
% %lower connectors
% for a = 1:m
% 	SXYnp(a + (l1-1)*Lx, a + 1 + l1*Lx) = -1/4;
% 	SXYnp(a + 1 + l1*Lx, a + (l1-1)*Lx) = -1/4;
% endfor
% for a = 1:(Lx-m-1)
% 	SXYnp(a + m + (l1-1)*Lx, a + m + 2 + l1*Lx) = -1/4;
% 	SXYnp(a + m + 2 + l1*Lx, a + m + (l1-1)*Lx) = -1/4;
% endfor
% for a = 1:m-1
% 	SXYnp(a + 1 + (l1-1)*Lx, a + l1*Lx) = 1/4;
% 	SXYnp(a + l1*Lx, a + 1 + (l1-1)*Lx) = 1/4;
% endfor
% for a = 1:(Lx-m)
% 	SXYnp(a + m + (l1-1)*Lx, a + m + l1*Lx) = 1/4;
% 	SXYnp(a + m + l1*Lx, a + m + (l1-1)*Lx) = 1/4;
% endfor

% %upper connectors
% for a = 1:(m-1)
% 	SXYnp(a + Lx*(l1+l2) + l2 - (Lx + 1), a + 1 + Lx*(l1+l2) + l2) = -1/4;
% 	SXYnp(a + 1 + Lx*(l1+l2) + l2, a + Lx*(l1+l2) + l2 - (Lx + 1)) = -1/4;
% endfor
% for a = 1:((Lx + 1) - (m+1))
% 	SXYnp(a + m + Lx*(l1 + l2) + l2 - (Lx + 1), a + m + Lx*(l1 + l2) + l2) = -1/4;
% 	SXYnp(a + m + Lx*(l1 + l2) + l2, a + m + Lx*(l1 + l2) + l2 - (Lx + 1)) = -1/4;
% endfor
% for a = 1:m
% 	SXYnp(a + 1 + Lx*(l1 + l2) + l2 - (Lx + 1), a + Lx*(l1 + l2) + l2) = 1/4;
% 	SXYnp(a + Lx*(l1 + l2) + l2, a + 1 + Lx*(l1+l2) + l2 - (Lx + 1)) = 1/4;
% endfor
% for a = 1:(Lx - m - 1)
% 	SXYnp(a + m + 2 + Lx*(l1 + l2) + l2 - (Lx + 1), a + m + Lx*(l1 + l2) + l2) = 1/4;
% 	SXYnp(a + m + Lx*(l1 + l2) + l2, a + m + 2 + Lx*(l1 + l2) + l2 - (Lx + 1)) = 1/4;
% endfor

% %%%%Periodic Boundary Condition%%%%%
% SXYp = SXYnp;
% for b = 1:l1
% 	SXYp(b * Lx, b * Lx + 1) = -1/4;
% 	SXYp(b * Lx + 1, b * Lx) = -1/4;
% endfor
% for b = 1:l2
% 	SXYp(b*(Lx + 1) + l1 * Lx, b*(Lx+1) + l1 * Lx + 1) = -1/4;
% 	SXYp(b*(Lx + 1) + l1 * Lx + 1, b*(Lx+1) + l1 * Lx) = -1/4;
% endfor
% for b = 1:(Ly - (l1 + l2) - 1)
% 	SXYp(b*Lx + (l1 + l2)*Lx + l2, b*Lx + (l1 + l2)*Lx + l2 + 1) = -1/4;
% 	SXYp(b*Lx + (l1 + l2)*Lx + l2 + 1, b*Lx + (l1 + l2)*Lx + l2) = -1/4;
% endfor
% for b = 1:(l1 - 1)
% 	SXYp((b+1)*Lx, (b-1)*Lx + 1) = 1/4;
% 	SXYp((b-1)*Lx + 1, (b+1)*Lx) = 1/4;
% endfor

% SXYp((l1+1)*Lx + 1, (l1 - 1)*Lx + 1) = 1/4;
% SXYp((l1 - 1)*Lx + 1, (l1+1)*Lx + 1) = 1/4;

% for b = 1:(l2-1)
% 	SXYp((b+1)*(Lx + 1) + l1*Lx, (Lx + 1)*(b-1) + l1*Lx + 1) = 1/4;
% 	SXYp((Lx + 1)*(b-1) + l1*Lx + 1, (b+1)*(Lx + 1) + l1*Lx) = 1/4;
% endfor

% SXYp((l1 + l2 + 1)*Lx + l2, (l1 + l2)*Lx + l2 - (Lx + 1)) = 1/4;
% SXYp((l1 + l2)*Lx + l2 - (Lx + 1), (l1 + l2 + 1)*Lx + l2) = 1/4;

% for b = 1:(Ly - (l1 + l2) -1)
% 	SXYp((b+1)*Lx + (l1 + l2)*Lx + l2, (b-1)*Lx + (l1 + l2)*Lx + l2 - 1) = 1/4;
% 	SXYp((b-1)*Lx + (l1 + l2)*Lx + l2 - 1, (b+1)*Lx + (l1 + l2)*Lx + l2) = 1/4;
% endfor

% for a = 1:(Lx - 1)
% 	SXYp((Ly - 1)*Lx + l2 + a, a+1) = -1/4;
% 	SXYp(a+1, (Ly - 1)*Lx + l2 + a) = -1/4;

% 	SXYp((Ly-1)*Lx + l2 + a + 1, a) = 1/4;
% 	SXYp(a, (Ly-1)*Lx + l2 + a + 1) = 1/4;
% endfor

% %Corners
% SXYp(Lx, Lx*Ly + l2 +1 - Lx) = 1/4;
% SXYp(Lx*Ly + l2 +1 - Lx, Lx) = 1/4;

% SXYp(Lx*Ly + l2, 1) = -1/4;
% SXYp(1, Lx*Ly + l2) = -1/4;
