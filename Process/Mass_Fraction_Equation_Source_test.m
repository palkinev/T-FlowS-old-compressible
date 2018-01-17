% octave script to test Mass_Fraction_Equation_Source
% to use this function you need to edit
% 1) Mass_Fraction_Equation_Source.f90 to make output
% 2) inivar to set specific time and y (xc(c), 0.03125, 0.125)
% 3) Processor.f90 so Mass_Fraction_Equation_Source is called
%    just after inivar.f90

% check on inflow and outflow!

rho0 = 5;
rho1 = 1;
k = 2.0;
w = 2.0;
vis = 0.001;
uf = 0.5;

Nx=32;
dx = 2/Nx; % dx = dy, dz = 2
dV = 2*dx^2;
dl = dV^(2/3); % ((dV)^(1/3))^2
x = -1 + dx/2 : dx : 1 - dx/2;
y = 0.03125;
t = 0.125;

xt = x - uf*t;
yt = y - uf*t;

z_ex = (1.0 + sin (k * pi * xt).*sin (k * pi * yt).*cos (w * pi * t))./ ...
     (1.0 + rho0/rho1+(1.0-rho0/rho1) .* sin (k * pi * xt) .* sin (k * pi * yt) .* cos (w * pi * t));

rho_ex = 1.0 ./(z_ex/rho1 + (1.0-z_ex)/rho0);

u_ex = -w/k/4.0 * cos (k * pi * xt) .* sin (k * pi * yt) .* sin (w * pi * t) .* (rho1-rho0) ./ rho_ex;

% derevatives
tilde_z_x = zeros(1, size(rho_ex, 2)); % tilde_z_x
for i = 1:Nx
  switch(i)
    case(1)
      tilde_z_x(i) = ( 1/2*z_ex(i+1) - 1/2*z_ex(Nx)  ) / dx;
    case(Nx)
      tilde_z_x(i) = ( 1/2*z_ex(1)   - 1/2*z_ex(i-1) ) / dx;
    otherwise
      tilde_z_x(i) = ( 1/2*z_ex(i+1) - 1/2*z_ex(i-1) ) / dx;
    endswitch
endfor

% derevatives x 2
tilde_z_xx = zeros(1, size(rho_ex, 2)); % tilde_z_xx
for i = 1:Nx
  switch(i)
    case(1)
      tilde_z_xx(i) = ( 1/2*tilde_z_x(i+1) - 1/2*tilde_z_x(Nx)  ) / dx;
    case(Nx)
      tilde_z_xx(i) = ( 1/2*tilde_z_x(1)   - 1/2*tilde_z_x(i-1) ) / dx;
    otherwise
      tilde_z_xx(i) = ( 1/2*tilde_z_x(i+1) - 1/2*tilde_z_x(i-1) ) / dx;
    endswitch
endfor

lap_ex = z_ex - dl * (tilde_z_xx.^2)/24;

%comparing with T-flowS
fid = fopen ("/home/l_palkin_e/eclipse/compressible/eclipse-project/Test_cases/Problem3/32/lap.dat", "r");
A = transpose(fscanf(fid,"%e",[18 Inf]));
fclose(fid);



x_f = A(:,1)';
y_f = A(:,2)';
z_f = A(:,3)';
u_f = A(:,4)';
u_ex_f = A(:,5)';
v_f = A(:,6)';
v_ex_f = A(:,7)';
w_f = A(:,8)';
z_f = A(:,9)';
z_ex_f = A(:,10)';
rho_f = A(:,11)';
rho_ex_f = A(:,12)';
lap_f = A(:,13)';
tilde_z_x_f = A(:,14)';
tilde_z_xx_f = A(:,15)';
tilde_z_yy_f = A(:,16)';
tilde_z_zz_f = A(:,17)';
dl_f = A(:,18)';

figure(1);
subplot(3,3,1)
plot (x, abs(z_f-z_ex),'ro-', 'LineWidth', 2,'MarkerSize',5); title ("|Z_{ex}-Z_{file}|"); xlabel('x, y=\Delta x/2'); ylabel('tolerance'); grid on;
subplot(3,3,2)
plot (x, abs(u_f-u_ex),'ro-', 'LineWidth', 2,'MarkerSize',5); title ("|u_{ex}-u_{file}|"); xlabel('x, y=\Delta x/2'); ylabel('tolerance'); grid on;
subplot(3,3,3)
plot (x, abs(rho_f-rho_ex),'ro-', 'LineWidth', 2,'MarkerSize',5); title ("|\\rho_{ex}-\\rho_{file}|"); xlabel('x, y=\Delta x/2'); ylabel('tolerance'); grid on;
subplot(3,3,4)
plot (x, abs(dl_f-dl),'ro-', 'LineWidth', 2,'MarkerSize',5); title ("|dl_{ex}-dl_{file}|"); xlabel('x, y=\Delta x/2'); ylabel('tolerance'); grid on;
subplot(3,3,5)
plot (x, abs(tilde_z_x_f-tilde_z_x),'ro-', 'LineWidth', 2,'MarkerSize',5); title ("|dz/dx_{ex}-dz/dx_{file}|"); xlabel('x, y=\Delta x/2'); ylabel('tolerance'); grid on;
subplot(3,3,6)
plot (x, abs(tilde_z_xx_f-tilde_z_xx),'ro-', 'LineWidth', 2,'MarkerSize',5); title ("|d2z/dx2_{ex}-d2z/dx2_{file}|"); xlabel('x, y=\Delta x/2'); ylabel('tolerance'); grid on;
subplot(3,3,7)
plot(x, abs(tilde_z_yy_f),'ro-', 'LineWidth', 2,'MarkerSize',5); title ("|d2z/dy2_{ex}-d2z/dy2_{file}|"); xlabel('x, y=\Delta x/2'); ylabel('tolerance'); grid on;
subplot(3,3,8)
plot(x, abs(tilde_z_zz_f),'ro-', 'LineWidth', 2,'MarkerSize',5); title ("|d2z/dz2_{ex}-d2z/dz2_{file}|"); xlabel('x, y=\Delta x/2'); ylabel('tolerance'); grid on;
subplot(3,3,9)
plot(x, abs(lap_f-lap_ex),'ro-', 'LineWidth', 2,'MarkerSize',5); title ("|new\\_Z_{ex}-new\\_Z_{file}|"); xlabel('x, y=\Delta x/2'); ylabel('tolerance'); grid on;