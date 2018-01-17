% octave script to test Calc_Sgs_Coefficients_Dynamic_Smagorinsky
% to use this function you need to edit
% 1) Calc_Sgs_Coefficients_Dynamic_Smagorinsky.f90 to make output
% 2) inivar to set specific time and y (xc(c), 0.03125, 0.125)
% 3) Processor.f90 so Calc_Sgs_Coefficients_Dynamic_Smagorinsky is called
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
tilde_u_x = zeros(1, size(rho_ex, 2)); % tilde_u_x
for i = 1:Nx
  switch(i)
    case(1)
      %tilde_u_x(i) = ( - 3/2*u_ex(i) + 2*u_ex(i+1) - 1/2*u_ex(i+2) ) / dx;
      tilde_u_x(i) = ( 1/2*u_ex(i+1) - 1/2*u_ex(Nx)  ) / dx;
    case(Nx)
      %tilde_u_x(i) = (  3/2*u_ex(i) - 2*u_ex(i-1) + 1/2*u_ex(i-2) ) / dx;
      tilde_u_x(i) = ( 1/2*u_ex(1)   - 1/2*u_ex(i-1) ) / dx;
    otherwise
      tilde_u_x(i) = ( 1/2*u_ex(i+1) - 1/2*u_ex(i-1) ) / dx;
    endswitch
endfor

tilde_sij_11 =       tilde_u_x;

pi_g = tilde_sij_11.^2.;

hat_volume           = zeros(1, size(rho_ex, 2)); % hat_volume
hat_rho              = zeros(1, size(rho_ex, 2)); % hat_rho
hat_rui_1            = zeros(1, size(rho_ex, 2)); % hat_rui_1
hat_r_tilde_sij_11   = zeros(1, size(rho_ex, 2)); % hat_r_tilde_sij_11
hat_ruiuj_11         = zeros(1, size(rho_ex, 2)); % hat_ruiuj_11
hat_p_term           = zeros(1, size(rho_ex, 2)); % hat_p_term
hat_mij_term_11      = zeros(1, size(rho_ex, 2)); % hat_mij_term_11


for i = 1:Nx
  switch(i)
    case(1)
      Nim1 = 1;
      Nip1 = 1;
      hat_volume        (i) = (    dV                                              *Nim1 +     dV                                        *5 +     dV                                              *Nip1);
      hat_rho           (i) = (    rho_ex                                     (Nx )*Nim1 +     rho_ex                                 (i)*5 +     rho_ex                                     (i+1)*Nip1);
      hat_rui_1         (i) = (    rho_ex(Nx )*u_ex                           (Nx )*Nim1 +     rho_ex(i)*u_ex                         (i)*5 +     rho_ex(i+1)*u_ex                           (i+1)*Nip1);
      hat_r_tilde_sij_11(i) = (    rho_ex(Nx )*tilde_sij_11                   (Nx )*Nim1 +     rho_ex(i)*tilde_sij_11                 (i)*5 +     rho_ex(i+1)*tilde_sij_11                   (i+1)*Nip1);
      hat_ruiuj_11      (i) = (    rho_ex(Nx )*u_ex(Nx )*u_ex                 (Nx )*Nim1 +     rho_ex(i)*u_ex(i)*u_ex                 (i)*5 +     rho_ex(i+1)*u_ex(i+1)*u_ex                 (i+1)*Nip1);
      hat_p_term        (i) = (  2*rho_ex(Nx )*dl*pi_g                        (Nx )*Nim1 +   2*rho_ex(i)*dl*pi_g                      (i)*5 +   2*rho_ex(i+1)*dl*pi_g                        (i+1)*Nip1);
      hat_mij_term_11   (i) = (4/3*rho_ex(Nx )*dl*pi_g(Nx )^(0.5)*tilde_sij_11(Nx )*Nim1 + 4/3*rho_ex(i)*dl*pi_g(i)^(0.5)*tilde_sij_11(i)*5 + 4/3*rho_ex(i+1)*dl*pi_g(i+1)^(0.5)*tilde_sij_11(i+1)*Nip1);
      %hat_rho           (i) = (                                                                rho_ex                                 (i)*5 +     rho_ex                                     (i+1)*Nip1);
    case(Nx)
      Nim1 = 1;
      Nip1 = 1;
      hat_volume        (i) = (    dV                                              *Nim1 +     dV                                        *5 +     dV                                              *Nip1);
      hat_rho           (i) = (    rho_ex                                     (i-1)*Nim1 +     rho_ex                                 (i)*5 +     rho_ex                                     (1  )*Nip1);
      hat_rui_1         (i) = (    rho_ex(i-1)*u_ex                           (i-1)*Nim1 +     rho_ex(i)*u_ex                         (i)*5 +     rho_ex(1  )*u_ex                           (1  )*Nip1);
      hat_r_tilde_sij_11(i) = (    rho_ex(i-1)*tilde_sij_11                   (i-1)*Nim1 +     rho_ex(i)*tilde_sij_11                 (i)*5 +     rho_ex(1  )*tilde_sij_11                   (1  )*Nip1);
      hat_ruiuj_11      (i) = (    rho_ex(i-1)*u_ex(i-1)*u_ex                 (i-1)*Nim1 +     rho_ex(i)*u_ex(i)*u_ex                 (i)*5 +     rho_ex(1  )*u_ex(1  )*u_ex                 (1  )*Nip1);
      hat_p_term        (i) = (  2*rho_ex(i-1)*dl*pi_g                        (i-1)*Nim1 +   2*rho_ex(i)*dl*pi_g                      (i)*5 +   2*rho_ex(1  )*dl*pi_g                        (1  )*Nip1);
      hat_mij_term_11   (i) = (4/3*rho_ex(i-1)*dl*pi_g(i-1)^(0.5)*tilde_sij_11(i-1)*Nim1 + 4/3*rho_ex(i)*dl*pi_g(i)^(0.5)*tilde_sij_11(i)*5 + 4/3*rho_ex(1  )*dl*pi_g(1  )^(0.5)*tilde_sij_11(1  )*Nip1);
    otherwise
      Nim1 = 1;
      Nip1 = 1;
      hat_volume        (i) = (    dV                                              *Nim1 +     dV                                        *5 +     dV                                              *Nip1);
      hat_rho           (i) = (    rho_ex                                     (i-1)*Nim1 +     rho_ex                                 (i)*5 +     rho_ex                                     (i+1)*Nip1);
      hat_rui_1         (i) = (    rho_ex(i-1)*u_ex                           (i-1)*Nim1 +     rho_ex(i)*u_ex                         (i)*5 +     rho_ex(i+1)*u_ex                           (i+1)*Nip1);
      hat_r_tilde_sij_11(i) = (    rho_ex(i-1)*tilde_sij_11                   (i-1)*Nim1 +     rho_ex(i)*tilde_sij_11                 (i)*5 +     rho_ex(i+1)*tilde_sij_11                   (i+1)*Nip1);
      hat_ruiuj_11      (i) = (    rho_ex(i-1)*u_ex(i-1)*u_ex                 (i-1)*Nim1 +     rho_ex(i)*u_ex(i)*u_ex                 (i)*5 +     rho_ex(i+1)*u_ex(i+1)*u_ex                 (i+1)*Nip1);
      hat_p_term        (i) = (  2*rho_ex(i-1)*dl*pi_g                        (i-1)*Nim1 +   2*rho_ex(i)*dl*pi_g                      (i)*5 +   2*rho_ex(i+1)*dl*pi_g                        (i+1)*Nip1);
      hat_mij_term_11   (i) = (4/3*rho_ex(i-1)*dl*pi_g(i-1)^(0.5)*tilde_sij_11(i-1)*Nim1 + 4/3*rho_ex(i)*dl*pi_g(i)^(0.5)*tilde_sij_11(i)*5 + 4/3*rho_ex(i+1)*dl*pi_g(i+1)^(0.5)*tilde_sij_11(i+1)*Nip1);
    endswitch
%      hat_volume        (i) = hat_volume        (i)                ;
      hat_rho           (i) = hat_rho           (i) / (5+Nim1+Nip1);
      hat_rui_1         (i) = hat_rui_1         (i) / (5+Nim1+Nip1);
      hat_r_tilde_sij_11(i) = hat_r_tilde_sij_11(i) / (5+Nim1+Nip1);
      hat_ruiuj_11      (i) = hat_ruiuj_11      (i) / (5+Nim1+Nip1);
      hat_p_term        (i) = hat_p_term        (i) / (5+Nim1+Nip1);
      hat_mij_term_11   (i) = hat_mij_term_11   (i) / (5+Nim1+Nip1);
endfor

pi_t = hat_r_tilde_sij_11.^2 ./ hat_rho.^2;

lij_11 = hat_ruiuj_11 - hat_rui_1 .* hat_rui_1 ./ hat_rho;
lij_kk = lij_11;

hat_delta = hat_volume.^(1./3.);

mij_11 = hat_mij_term_11      - 2. * hat_delta.^2. .* pi_t.^(0.5) .* ( hat_r_tilde_sij_11*(2/3 ) );
mij_22 = hat_mij_term_11/(-2) - 2. * hat_delta.^2. .* pi_t.^(0.5) .* ( hat_r_tilde_sij_11*(-1/3) );
mij_33 = hat_mij_term_11/(-2) - 2. * hat_delta.^2. .* pi_t.^(0.5) .* ( hat_r_tilde_sij_11*(-1/3) );
mij_kk = mij_11 + mij_22 + mij_33;

mijmij = mij_11.^2 + mij_22.^2 + mij_33.^2;

mijldij = mij_11 .* lij_11 - mij_kk.*lij_kk/3;

c_i = lij_kk ./ ( hat_p_term - 2. * hat_delta.^2. .* hat_rho .* pi_t );

c_r = mijldij ./ ( mijmij );

figure(1);
subplot(1,3,1)
plot (x, rho_ex.*tilde_sij_11,'ro-', x, hat_r_tilde_sij_11,'bx-')
subplot(1,3,2)
plot (x, rho_ex,'ro-', x, hat_rho,'bx-')
subplot(1,3,3)
plot (x, u_ex.*rho_ex,'ro-', x, hat_rui_1,'bx-')


%comparing with T-flowS
fid = fopen ("/home/l_palkin_e/eclipse/compressible/eclipse-project/Test_cases/Problem3/32/hat.dat", "r");
A = transpose(fscanf(fid,"%e",[66 Inf]));
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
hat_volume_f = A(:,13)';
hat_rho_f = A(:,14)';
hat_rui_1_f = A(:,15)';

hat_r_tilde_sij_11_f = A(:,18)';
hat_r_tilde_sij_22_f = A(:,19)';
hat_r_tilde_sij_33_f = A(:,20)';
hat_r_tilde_sij_12_f = A(:,21)';
hat_r_tilde_sij_13_f = A(:,22)';
hat_r_tilde_sij_23_f = A(:,23)';

hat_ruiuj_11_f = A(:,24)';

hat_p_term_f = A(:,30)';
hat_mij_term_11_f = A(:,31)';
pi_t_f = A(:,37)';
lij_11_f = A(:,38)';
lij_22_f = A(:,39)';
lij_33_f = A(:,40)';
lij_12_f = A(:,41)';
lij_13_f = A(:,42)';
lij_23_f = A(:,43)';
hat_delta_f = A(:,44)';

mij_11_f = A(:,45)';
mij_22_f = A(:,46)';
mij_33_f = A(:,47)';
mij_12_f = A(:,48)';
mij_13_f = A(:,49)';
mij_23_f = A(:,50)';

mijmij_f = A(:,51)';
mijldij_f = A(:,52)';


c_i_f = A(:,53)';
c_r_f = A(:,54)';
tilde_u_x_f = A(:,55)';
tilde_u_y_f = A(:,56)';
tilde_u_z_f = A(:,57)';
tilde_v_x_f = A(:,58)';
tilde_v_y_f = A(:,59)';
tilde_v_z_f = A(:,60)';
tilde_w_x_f = A(:,61)';
tilde_w_y_f = A(:,62)';
tilde_w_z_f = A(:,63)';
vol_f = A(:,64)';
pi_g_f = A(:,65)';
bar_delta_f = A(:,66)';

figure(2);
subplot(5,4,1)
plot(x, abs(u_ex-u_f),'ro-', 'LineWidth', 2,'MarkerSize',5); title ("|u_{ex}-u_{file}|"); xlabel('x, y=\Delta x/2'); ylabel('tolerance'); grid on;
subplot(5,4,2)
plot(x, abs(hat_rho-hat_rho_f),'ro-', 'LineWidth', 2,'MarkerSize',5); title ("|\\rho_{ex}^f-\\rho_{file}^f|"); xlabel('x, y=\Delta x/2'); ylabel('tolerance'); grid on;
subplot(5,4,3)
plot(x, abs(hat_volume-hat_volume_f),'ro-', 'LineWidth', 2,'MarkerSize',5); title ("|V_{ex}^f-V_{file}^f|"); xlabel('x, y=\Delta x/2'); ylabel('tolerance'); grid on;
subplot(5,4,4)
plot(x, abs(rho_ex-rho_f),'ro-', 'LineWidth', 2,'MarkerSize',5); title ("|\\rho_{ex}-\\rho_{file}|"); xlabel('x, y=\Delta x/2'); ylabel('tolerance'); grid on;
subplot(5,4,5)
plot(x, abs(tilde_u_x-tilde_u_x_f),'ro-', 'LineWidth', 2,'MarkerSize',5); title ("|du/dx_{ex}-du/dx_{file}|"); xlabel('x, y=\Delta x/2'); ylabel('tolerance'); grid on;
subplot(5,4,6)
plot(x, abs(vol_f-(dx^2)*2),'ro-', 'LineWidth', 2,'MarkerSize',5); title ("|vol_{ex}-vol_{file}|"); xlabel('x, y=\Delta x/2'); ylabel('tolerance'); grid on;
subplot(5,4,7)
plot(x, abs(hat_rui_1-hat_rui_1_f),'ro-', 'LineWidth', 2,'MarkerSize',5); title ("|rui\\_1_{ex}^f-rui\\_1_{file}^f|"); xlabel('x, y=\Delta x/2'); ylabel('tolerance'); grid on;
subplot(5,4,8)
plot(x, abs(hat_r_tilde_sij_11-hat_r_tilde_sij_11_f),'ro-', 'LineWidth', 2,'MarkerSize',5); title ("|r\\_sij\\_11_{ex}^f-r\\_sij\\_11_{file}^f|"); xlabel('x, y=\Delta x/2'); ylabel('tolerance'); grid on;
subplot(5,4,9)
plot(x, abs(hat_ruiuj_11-hat_ruiuj_11_f),'ro-', 'LineWidth', 2,'MarkerSize',5); title ("|ruiuj\\_11_{ex}^f-ruiuj\\_11_{file}^f|"); xlabel('x, y=\Delta x/2'); ylabel('tolerance'); grid on;
subplot(5,4,10)
plot(x, abs(pi_g-pi_g_f),'ro-', 'LineWidth', 2,'MarkerSize',5); title ("|\\Pi_G_{ex}-\\Pi_G_{file}|"); xlabel('x, y=\Delta x/2'); ylabel('tolerance'); grid on;
subplot(5,4,11)
plot(x, abs(pi_t-pi_t_f),'ro-', 'LineWidth', 2,'MarkerSize',5); title ("|\\Pi_T_{ex}^f-\\Pi_T_{file}^f|"); xlabel('x, y=\Delta x/2'); ylabel('tolerance'); grid on;
subplot(5,4,12)
plot(x, abs(hat_p_term-hat_p_term_f),'ro-', 'LineWidth', 2,'MarkerSize',5); title ("|p\\_term_{ex}^f-p\\_term_{file}^f|"); xlabel('x, y=\Delta x/2'); ylabel('tolerance'); grid on;
subplot(5,4,13)
plot(x, abs(mij_11-mij_11_f),'ro-', 'LineWidth', 2,'MarkerSize',5); title ("|mij\\_11_{ex}^f-mij\\_11_{file}^f|"); xlabel('x, y=\Delta x/2'); ylabel('tolerance'); grid on;
subplot(5,4,14)
plot(x, abs(lij_11-lij_11_f),'ro-', 'LineWidth', 2,'MarkerSize',5); title ("|lij\\_11_{ex}^f-lij\\_11_{file}^f|"); xlabel('x, y=\Delta x/2'); ylabel('tolerance'); grid on;
subplot(5,4,15)
plot(x, abs(mijmij-mijmij_f),'ro-', 'LineWidth', 2,'MarkerSize',5); title ("|mijmij_{ex}^f-mijmij_{file}^f|"); xlabel('x, y=\Delta x/2'); ylabel('tolerance'); grid on;
subplot(5,4,16)
plot(x, abs(mijldij-mijldij_f),'ro-', 'LineWidth', 2,'MarkerSize',5); title ("|mijldij_{ex}^f-mijmij_{file}^f|"); xlabel('x, y=\Delta x/2'); ylabel('tolerance'); grid on;
subplot(5,4,19)
plot(x, abs(c_i-c_i_f),'ro-', 'LineWidth', 2,'MarkerSize',5); title ("|c\\_i_{ex}^f-c\\_i_{file}^f|"); xlabel('x, y=\Delta x/2'); ylabel('tolerance'); grid on;
subplot(5,4,20)
plot(x, abs(c_r-c_r_f),'ro-', 'LineWidth', 2,'MarkerSize',5); title ("|c\\_r_{ex}^f-c\\_r_{file}^f|"); xlabel('x, y=\Delta x/2'); ylabel('tolerance'); grid on;
%plot(x, abs(vol_f-(dx^2)*2),'ro-', 'LineWidth', 2,'MarkerSize',5); title ("|vol_{ex}-vol_{file}|"); xlabel('x, y=\Delta x/2'); ylabel('tolerance'); grid on;

