m = 5000;
k = 3000;
wn = sqrt(k/m);
zeta = 0.02;
c = 2*zeta*wn*m;

K = [2*k -k 0 0 0;  %stifness matrix
	-k 2*k -k 0 0;
	0 -k 2*k -k 0;
	0 0 -k 2*k -k;
	0 0 0 -k k];

C = [2*c -c 0 0 0; %damping matrix
	-c 2*c -c 0 0;
	0 -c 2*c -c 0;
	0 0 -c 2*c -c;
	0 0 0 -c c];

M = eye(5)*m; %initialising mass matrix


[phi, wsq] = eig(K, M); %eigenvalue function 


modal_damping = phi'*C*phi;  
modal_mass = phi'*M*phi;
modal_stiffness = phi'*K*phi;
