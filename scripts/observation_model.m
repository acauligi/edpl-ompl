K = sym('K', [3,3], 'real');  % intrinsics matrix
lk = sym('lk', [3,1], 'real');    % landmark position in world frame

% body to camera transform matrix 
Rbc  = sym('Rbc', [3,3], 'real');
pbc  = sym('pbc', [3,1], 'real');
Tbc  = [Rbc, pbc; zeros(1,3), 1];

x = sym('z', [12,1], 'real'); % flat output consisting of 
                              % (sigma, sigma_dot, sigma_ddot)
                              % where sigma = (x,y,z,yaw)
t = x(9:11); % acceleration component of state
psi = x(3);
zb = t/norm(t);

%% Method 1
xc = [cos(psi);sin(psi);0];
yb = cross(zb,xc);
xb = cross(yb,zb);
Rwb = [xb, yb, zb];     % body to world transformation matrix
Rbw = Rwb';
Hx1 = simplify(jacobian(Rbw(:,1),x), 'Steps', 100);
Hy1 = simplify(jacobian(Rbw(:,2),x), 'Steps', 100);
Hz1 = simplify(jacobian(Rbw(:,3),x), 'Steps', 100);

% fid = fopen('sensor_jacobian_Hx1.txt', 'w');
% fprintf(fid, char(Hx1)); 
% fclose(fid);
% fid = fopen('sensor_jacobian_Hy1.txt', 'w');
% fprintf(fid, char(Hy1)); 
% fclose(fid);
% fid = fopen('sensor_jacobian_Hz1.txt', 'w');
% fprintf(fid, char(Hz1)); 
% fclose(fid);

%% Method 2 
yc = [-sin(psi); cos(psi); 0];
xb = cross(yc,zb);
yb = cross(zb,xb);
Rwb = [xb, yb, zb];     % body to world transformation matrix
Rbw = Rwb';
Hx2 = simplify(jacobian(Rbw(:,1),x), 'Steps', 100);
Hy2 = simplify(jacobian(Rbw(:,2),x), 'Steps', 100);
Hz2 = simplify(jacobian(Rbw(:,3),x), 'Steps', 100);

fid = fopen('sensor_jacobian_Hx2.txt', 'w');
fprintf(fid, char(Hx2)); 
fclose(fid);
fid = fopen('sensor_jacobian_Hy2.txt', 'w');
fprintf(fid, char(Hy2)); 
fclose(fid);
fid = fopen('sensor_jacobian_Hz2.txt', 'w');
fprintf(fid, char(Hz2)); 
fclose(fid);
