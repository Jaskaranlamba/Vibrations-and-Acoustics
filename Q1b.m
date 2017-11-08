m = 5000;
k = 3000;
wn = sqrt(k/m);
si = 0.02;
c=2*si*wn*m;
dof = 5;


ga = 0.25; %gamma
be = 0.25; %beta
load elcentro.mat

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


%[phi, lambda] = eig(K, M); %eigenvalue function 
[phi,lamda]=eig(K,M);





modal_damping = phi'*C*phi;  
modal_mass = phi'*M*phi;
modal_stiffness = phi'*K*phi;

%%initial calculation
u = zeros(dof, length(t));
du = zeros(dof, length(t));
ddu = zeros(dof, length(t));
response = zeros(dof, length(t));
Po = [p;zeros(dof-1, length(t))];


u(:,1) =0;
du(:,1)=0;
response(:, 1) = phi*du(:, 1);
ddu(:, 1)=modal_mass\(p(:, 1)-modal_stiffness*u(:, 1)-modal_damping*du(:, 1));
dt= t(3) - t(2);
a = modal_mass/(be*dt) + (ga/be)*modal_damping;
b = modal_mass/(2*be) + (modal_damping*dt)*((ga/(2*be))-1);
kc = modal_stiffness +(ga/(be*dt))*modal_damping + (1/(be*dt^2))*modal_mass;


for j =1:1:length(t)-1
    P = phi*Po;
    dP = P(:,j+1)-P(:,j);
    dPo = dP+a*du(:,j)+b*ddu(:,j);
    inu = kc\dPo;
    indu = (inu*ga)/(be*dt)-(du(:,j)*ga)/be+dt*ddu(:,j)*(1-(ga/(2*be)));
    inddu = (1/(be*dt^2))*inu-(1/(be*dt))*du(:,j)-(1/(2*be))*ddu(:,j);
    u(:,j+1)= du(:,j)+inu;
    du(:,j+1)=du(:,j)+indu;
    ddu(:,j+1)=ddu(:,j)+inddu;
    response(:,j+1)=m*phi*u(:,j+1);
    
end


plot(t,response)
xlabel('t - time');
ylabel('u - response');
legend('m1', 'm2', 'm3', 'm4', 'm5');