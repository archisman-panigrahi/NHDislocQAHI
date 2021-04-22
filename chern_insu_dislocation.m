%%%% This is the main file %%%%
t = 1;
% For non-Hermitian matrix A, [b,c,d] = eig(A) gives
% The right eigenkets in columns of b
% The eigenvalues in the diagonal entries of c
% the kets corresponding to the left eigenbras in the columns of d
% For example, d(:,1)' * A = c(1,1) * d(:,1)'

%% Checks if a matrix is Hermitian
%% Prints 0 if it is Hermitian, and a non-zero number otherwise
%% Also prints the entries m,n if A(m,n) != A(n,m)^*
function Hermitian_check(a)
	b = size(a)
	c = 0;
	for xi = 1:b(1)
		for yi = 1:b(2)
			if (a(xi,yi) - conj(a(yi,xi))) == 0
				c += 0;
			else
				c += 1;
				xi,yi
			endif
		endfor
	endfor
	c
endfunction

%%%% saves LDoS to a given filename %%%%
function plotState(local_density, Lx, Ly, l1, l2, m, name)
	
	%lower left
	for a = 1:m
		for b = 1:l1
			density(a,b) = local_density(a + (b-1)*Lx);
		endfor
	endfor
	%lower gap
	for b = 1:l1
		density(m+1,b) = 0;
	endfor
	%lower right
	for a = (m+1):Lx
		for b = 1:l1
			density(a+1,b) = local_density(a + (b-1)*Lx);
		endfor
	endfor
	%middle part
	for a = 1:(Lx+1)
		for b = (l1 + 1):(l1 + l2)
			density(a,b) = local_density(l1 * Lx + a + (b-l1-1)*(Lx + 1));
		endfor
	endfor
	%upper left
	for a = 1:m
		for b = (l1 + l2 + 1):Ly
			density(a,b) = local_density(l1 * Lx + l2 * (Lx + 1) + a + (b-l1-l2-1)*Lx);
		endfor
	endfor
	%upper gap
	for b = (l1 + l2 + 1):Ly
		density(m+1,b) = 0;
	endfor
	%upper right
	for a = (m+1):Lx
		for b = (l1 + l2 + 1):Ly
			density(a+1,b) = local_density(l1 * Lx + l2 * (Lx + 1) + a + (b-l1-l2-1)*Lx);
		endfor
	endfor
	saveThisMatrixToFile(density, name);
endfunction


function NHenergySpectra(Lx,Ly,l1,l2,m,m_0,hx,hy,hz,h0)
	%This is for PBC
	close all;
	tic;
	angmom;
	BuildingBlocksDoubleDislocation;
	n_states = n_atoms * 2;
	h_NH = j * kron(M2D, hx * sigma_x + hy * sigma_y + hz * sigma_z + h0 * eye(2));
	h_p = kron(CXp + CYp - m_0*M2D, sigma_z) + kron(SXp, sigma_x) + kron(SYp, sigma_y) + h_NH;

	[states_p, energies_p, left_p] = eig(h_p);

	%%%%Arrange energies in ascending order of their real part
	energies_p = diag(energies_p);

	[~,I_p] = sort(real(energies_p));

	energies_p = energies_p(I_p);

	states_p = states_p(:,I_p);
	left_p = left_p(:,I_p);

	for a = 1:n_atoms
		% Each site can host two states. Here we are taking 
		% |<site_state_1| psi_1>|^2 + |<site_state_2| psi_1>|^2 + |<site_state_1| psi_2>|^2 + |<site_state_2| psi_2>|^2
		that_state_p(a) = (abs(states_p(2*a - 1, n_atoms))).^2 + (abs(states_p(2*a - 1, n_atoms + 1))).^2 + (abs(states_p(2*a, n_atoms))).^2 + (abs(states_p(2*a, n_atoms + 1))).^2;
		that_state_left_p(a) = (abs(left_p(2*a - 1, n_atoms))).^2 + (abs(left_p(2*a - 1, n_atoms + 1))).^2 + (abs(left_p(2*a, n_atoms))).^2 + (abs(left_p(2*a, n_atoms + 1))).^2;
		that_state_biorthogonal_p(a) = abs(left_p(2*a - 1, n_atoms) * states_p(2*a - 1, n_atoms)) + abs(left_p(2*a - 1, n_atoms + 1) * states_p(2*a - 1, n_atoms + 1)) + abs(left_p(2*a, n_atoms) * states_p(2*a, n_atoms)) + abs(left_p(2*a, n_atoms + 1) * states_p(2*a, n_atoms + 1));
	endfor

	%%%%% Saves the densities to .dat files
	name2 = strcat('Lx_',num2str(Lx),'Ly_',num2str(Ly),'l1_',num2str(l1),'l2_',num2str(l2),'m_',num2str(m))

	cd saved_plots/NHdoubleDislocationChern/states
	dirname1 = strcat('m',num2str(m_0));
	dirname2 = strcat('hx',num2str(hx),'hy',num2str(hy),'hz',num2str(hz));
	mkdir(num2str(dirname1));
	cd(num2str(dirname1));
	mkdir(num2str(dirname2));
	cd ../../../../

	plotState(that_state_p,Lx,Ly,l1,l2,m,strcat('PBC_state',name2));
	plotState(that_state_left_p,Lx,Ly,l1,l2,m,strcat('PBC_state',name2,'left'));
	plotState(that_state_biorthogonal_p,Lx,Ly,l1,l2,m,strcat('PBC_state',name2,'biorthogonal'));

	system(strcat('mv',"\t", 'PBC_state', name2, '.dat', "\t", 'saved_plots/NHdoubleDislocationChern/states/', dirname1, '/', dirname2));
	system(strcat('mv',"\t", 'PBC_state', name2, 'left', '.dat', "\t", 'saved_plots/NHdoubleDislocationChern/states/', dirname1, '/', dirname2));
	system(strcat('mv',"\t", 'PBC_state', name2, 'biorthogonal', '.dat', "\t", 'saved_plots/NHdoubleDislocationChern/states/', dirname1, '/', dirname2));


	energies_p_export = [((real(energies_p))')(:), ((imag(energies_p))')(:)];
	dlmwrite(strcat('PBC_energies',name2,'.dat'), energies_p_export, 'delimiter', ' ');
	figure()
	scatter(real(energies_p),imag(energies_p))
	cd saved_plots/NHdoubleDislocationChern/energies
	mkdir(num2str(dirname1));
	cd(num2str(dirname1));
	mkdir(num2str(dirname2));
	cd ../../../..
	system(strcat('mv',"\t", 'PBC_energies', name2, '.dat', "\t", 'saved_plots/NHdoubleDislocationChern/energies/', dirname1, '/',dirname2));

	toc;
endfunction

function NHenergySpectraOBC(Lx,Ly,l1,l2,m,m_0,hx,hy,hz,h0)
	close all;
	tic;
	angmom;
	BuildingBlocksDoubleDislocation;
	n_states = n_atoms * 2;
	h_NH = j * kron(M2D, hx * sigma_x + hy * sigma_y + hz * sigma_z + h0 * eye(2));
	h_np = kron(CXnp + CYnp - m_0*M2D, sigma_z) + kron(SXnp, sigma_x) + kron(SYnp, sigma_y) + h_NH;

	[states_np, energies_np, left_np] = eig(h_np);

	%%%%Arrange energies in ascending order of their real part
	energies_np = diag(energies_np);

	[~,I_np] = sort(real(energies_np));

	energies_np = energies_np(I_np);

	states_np = states_np(:,I_np);
	left_np =	left_np(:,I_np);


	for a = 1:n_atoms
		% Each site can host two states. Here we are taking 
		% |<site_state_1| psi_1>|^2 + |<site_state_2| psi_1>|^2 + |<site_state_1| psi_2>|^2 + |<site_state_2| psi_2>|^2
		that_state_np(a) = (abs(states_np(2*a - 1, n_atoms))).^2 + (abs(states_np(2*a - 1, n_atoms + 1))).^2 + (abs(states_np(2*a, n_atoms))).^2 + (abs(states_np(2*a, n_atoms + 1))).^2;
		that_state_left_np(a) = (abs(left_np(2*a - 1, n_atoms))).^2 + (abs(left_np(2*a - 1, n_atoms + 1))).^2 + (abs(left_np(2*a, n_atoms))).^2 + (abs(left_np(2*a, n_atoms + 1))).^2;
		that_state_biorthogonal_np(a) = abs(left_np(2*a - 1, n_atoms) * states_np(2*a - 1, n_atoms)) + abs(left_np(2*a - 1, n_atoms + 1) * states_np(2*a - 1, n_atoms + 1)) + abs(left_np(2*a, n_atoms) * states_np(2*a, n_atoms)) + abs(left_np(2*a, n_atoms + 1) * states_np(2*a, n_atoms + 1));
	endfor

	%%%%% Saves the densities to .dat files
	name2 = strcat('Lx_',num2str(Lx),'Ly_',num2str(Ly),'l1_',num2str(l1),'l2_',num2str(l2),'m_',num2str(m))

	cd saved_plots/NHdoubleDislocationChern/states
	dirname1 = strcat('m',num2str(m_0));
	dirname2 = strcat('hx',num2str(hx),'hy',num2str(hy),'hz',num2str(hz));
	mkdir(num2str(dirname1));
	cd(num2str(dirname1));
	mkdir(num2str(dirname2));
	cd ../../../../
	
	plotState(that_state_np,Lx,Ly,l1,l2,m,strcat('OBC_state',name2));
	plotState(that_state_left_np,Lx,Ly,l1,l2,m,strcat('OBC_state',name2,'left'));
	plotState(that_state_biorthogonal_np,Lx,Ly,l1,l2,m,strcat('OBC_state',name2,'biorthogonal'));

	system(strcat('mv',"\t", 'OBC_state', name2, '.dat', "\t", 'saved_plots/NHdoubleDislocationChern/states/', dirname1, '/', dirname2))
	system(strcat('mv',"\t", 'OBC_state', name2, 'left', '.dat', "\t", 'saved_plots/NHdoubleDislocationChern/states/', dirname1, '/', dirname2))
	system(strcat('mv',"\t", 'OBC_state', name2, 'biorthogonal', '.dat', "\t", 'saved_plots/NHdoubleDislocationChern/states/', dirname1, '/', dirname2))

	energies_np_export = [((real(energies_np))')(:), ((imag(energies_np))')(:)];
	dlmwrite(strcat('OBC_energies',name2,'.dat'), energies_np_export, 'delimiter', ' ');

	figure()
	scatter(real(energies_np),imag(energies_np))
	axis tight;

	cd saved_plots/NHdoubleDislocationChern/energies
	mkdir(num2str(dirname1));
	cd(num2str(dirname1));
	mkdir(num2str(dirname2));
	cd ../../../..
	system(strcat('mv',"\t", 'OBC_energies', name2, '.dat', "\t", 'saved_plots/NHdoubleDislocationChern/energies/', dirname1, '/',dirname2))

	toc;
endfunction