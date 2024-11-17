#Constantes
G = 6.67408e-11;
Ma = 3.0e30;
Mb = 2.0e20;
##Mb = 2.0e28;
##Mb = 2.0e30; #Colision
Mc = 3.0e30;
AU = 1.5e11;
yearsec = 12*30*24.0 * 60 * 60; #Un año basicamente

gconst = G*[(Ma*Mb);(Mb*Mc); Mc*Ma];

#Condiciones iniciales

# a
r_a = [-2*AU, 0, 0];
v_a = [1e4, -1e4, 0];
r0_a = r_a;
v0_a = v_a;

#b
r_b = [0, 2*AU, 0];
v_b = [-1e4, -1e4, 0];
r0_b = r_b;
v0_b = v_b;

#c
r_c = [3*AU, 0, 0];
v_c = [-1e4, 1e4, 0];
r0_c = r_c;
v0_c = v_c;

t = 0.0;

dt = 0.02*yearsec; #Paso tiempo

alist = [];
blist = [];
clist = [];

valist = [];
vblist = [];
vclist = [];

#Comenzamos la simulacion
while t < 21 * yearsec
  r_ab = r_b - r_a;
  modr3_ab = norm(r_ab)^3;

  r_bc = r_c - r_b;
  modr3_bc = norm(r_bc)^3;

  r_ca = r_a - r_c;
  modr3_ca = norm(r_ca)^3;

  # a
  f_ab = (gconst(1)/modr3_ab)*r_ab;
  f_ac = -(gconst(3)/modr3_ca)*r_ca;

  v_a += ((f_ab + f_ac)*dt)/Ma;

  r_a += v_a*dt;

  alist(end +1, :) = r_a;
  valist(end +1, :) = v_a;

  # b
  f_ba = -(gconst(1)/modr3_ab)*r_ab;
  f_bc = (gconst(2)/modr3_bc)*r_bc;

  v_b += (f_ba + f_bc)*dt/Mb;

  r_b += v_b*dt;

  blist(end + 1, :) = r_b;
  vblist(end +1, :) = v_b;

  # c
  f_ca = (gconst(3)/modr3_ca)*r_ca;
  f_cb = -(gconst(2)/modr3_bc)*r_bc;

  v_c += (f_ca + f_cb)*dt/Mc;

  r_c += v_c*dt;

  clist(end + 1, :)  = r_c;
  vclist(end +1, :) = v_c;

  t = t + dt;

endwhile

disp("Data lista");

#Graficamos los resultados
figure(1);

plot3(alist(:, 1), alist(:, 2), alist(:, 3), 'b', blist(:, 1), blist(:, 2), blist(:, 3), 'r', clist(:, 1), clist(:, 2), clist(:, 3), 'g');
legend("Cuerpo A", "Cuerpo B", "Cuerpo C");

#Animacion para la simulacion de las trayectorias de los cuerpos
figure(2);
hold on;

#Definimos nuestras funciones que se encargaran de la animacion
function [h1, h2, h3] = Animate_2b(i, h1, h2, h3, r1_anim, r2_anim, r3_anim)
  delete(h1);
  delete(h2);
  delete(h3);

  plot3(r1_anim(1:i, 1), r1_anim(1:i, 2), r1_anim(1:i, 3), 'b');
  plot3(r2_anim(1:i, 1), r2_anim(1:i, 2), r2_anim(1:i, 3), 'r');
  plot3(r3_anim(1:i, 1), r3_anim(1:i, 2), r3_anim(1:i, 3), 'g');

  h1 = plot3(r1_anim(i, 1), r1_anim(i, 2), r1_anim(i, 3), 'bo', 'MarkerSize', 8);
  h2 = plot3(r2_anim(i, 1), r2_anim(i, 2), r2_anim(i, 3), 'r*', 'MarkerSize', 8);
  h3 = plot3(r3_anim(i, 1), r3_anim(i, 2), r3_anim(i, 3), 'go', 'MarkerSize', 8);
endfunction

function Animate(r1, r2, r3)
  #Creamos nuevos arrays para la animation
  r1_anim = r1(1:1:end, :);
  r2_anim = r2(1:1:end, :);
  r3_anim = r3(1:1:end, :);
  #Initial markers
  h1 = plot3(r1_anim(1, 1), r1_anim(1, 2), r1_anim(1, 3), 'bo', 'MarkerSize', 8, 'DisplayName', 'Cuerpo A');
  h2 = plot3(r2_anim(1, 1), r2_anim(1, 2), r2_anim(1,3), 'r*', 'MarkerSize', 8, 'DisplayName', 'Cuerpo B');
  h3 = plot3(r3_anim(1, 1), r3_anim(1, 2), r3_anim(1,3), 'go*', 'MarkerSize', 8, 'DisplayName', 'Cuerpo C');


  xlabel('Posicion X (m)');
  ylabel('Posicion Y (m)');
  zlabel('PosicionZ (m)');

  for i = 1:length(r1)
    [h1, h2, h3] = Animate_2b(i, h1, h2, h3, r1_anim, r2_anim, r3_anim);
    drawnow;
  endfor
endfunction

#Animamos las orbitas
Animate(alist, blist, clist);
title("Simulacion 3 Cuerpos");
xlim([-1.5e12, 1.5e12])
ylim([-1.5e12, 1.5e12])

#Pasamos ahora a interpolar las trayectorias usando spline

sp_a = spline(0:dt:(t-dt), alist, 0:7*24*60^2:(t-dt));
sp_b = spline(0:dt:(t-dt), blist, 0:7*24*60^2:(t-dt));
sp_c = spline(0:dt:(t-dt), clist, 0:7*24*60^2:(t-dt));

her_b = pchip(0:dt:(t-dt), blist', 0:7*24*60^2:(t-dt));

figure(3)
plot3(blist(:, 1), blist(:, 2), blist(:, 3),'r', sp_b(1, :), sp_b(2, :), sp_b(3, :), '-', her_b(1, :), her_b(2, :), her_b(3, :), 'go-');
legend("Cuerpo B", "Spline B", "Hermite B");
title("Trayectoria Cuerpo B Spline y Hermite vs Original");

figure(4)
plot3(alist(:, 1), alist(:, 2), alist(:, 3),'b', clist(:, 1), clist(:, 2), clist(:, 3),'g', sp_a(1, :), sp_a(2, :), sp_a(3, :), '-', sp_c(1, :), sp_c(2, :), sp_c(3, :), '-')
legend("Cuerpo A", "Cuerpo C", "Spline A", "Spline C");
title("Trayectoria Cuerpo A y C - Spline vs Original");

#Ahora usaremos tecnicas de diferenciacion e integracion numerica y ver que obtenemos

#Usando la formual del punto medio de tres puntos
dr_a = dr_b = dr_c = zeros(length(alist)-2, 3);

for i=2:length(alist)-1
  dr_a(i-1, :) = (alist(i+1, :) - alist(i-1,:))/(2*dt);
  dr_b(i-1, :) = (blist(i+1, :) - blist(i-1,:))/(2*dt);
  dr_c(i-1, :) = (clist(i+1, :) - clist(i-1,:))/(2*dt);
endfor

#Calculamos el error absoluto medio
err_a = sum(abs(valist(1:length(dr_a), :)-dr_a))/length(dr_a);
err_b = sum(abs(vblist(1:length(dr_a), :)-dr_b))/length(dr_b);
err_c = sum(abs(vclist(1:length(dr_a), :)-dr_c))/length(dr_c);

fprintf("Error Absoluto Medio Derivada\n")
fprintf("Error A %f | %f | %f \n", err_a);
fprintf("Error B %f | %f | %f \n", err_b);
fprintf("Error C %f | %f | %f \n", err_c);

# Con un h mas pequeño podriamos reducir drasticamente el error pero se eligio
#el dt que ya teniamos para asi facilitar la eleccion de puntos y la comparacion con nuestro vector v original
# Aun con el gran margen de error (dado que son numero e3, e4) podemos ver que
#los resultados con este metodo son similares al valor original. Ademas notamos que el error del
#cuerpo B es muchisimo mayor al de los otros dos cuerpos.

#Como usamos la formula de punto medio podemos notar que hacia el centro de los datos el error es mucho mas pequeño.
# Si usaramos la formula de punto extremo para evaluar los primeros y ultimos valores probablemente obtendriamos
#resultados con un menor error

#Hagamos el calculo para los primeros y ultimos 250 datos usando la formula de los extremos y evaluamos el error de nuevo
for i=1:250
  #Formula del extremo de tres puntos (extremo izquierdo)
  dr_a(i, :) = (-3*alist(i, :) + 4*alist(i+1, :) - alist(i+2,:))/(2*dt);
  dr_b(i, :) = (-3*blist(i, :) + 4*blist(i+1, :) - blist(i+2,:))/(2*dt);
  dr_c(i, :) = (-3*clist(i, :) + 4*clist(i+1, :) - clist(i+2,:))/(2*dt);

  #Formula del extremo de tres puntos (extremo derecho)
  dr_a(length(dr_a)+1-i, :) = (-3*alist(length(dr_a)+1-i, :) + 4*alist(length(dr_a)+1-(i+1), :) - alist(length(dr_a)+1-(i+2),:))/(-2*dt);
  dr_b(length(dr_b)+1-i, :) = (-3*blist(length(dr_b)+1-i, :) + 4*blist(length(dr_b)+1-(i+1), :) - blist(length(dr_b)+1-(i+2),:))/(-2*dt);
  dr_c(length(dr_c)+1-i, :) = (-3*clist(length(dr_c)+1-i, :) + 4*clist(length(dr_c)+1-(i+1), :) - clist(length(dr_c)+1-(i+2),:))/(-2*dt);
endfor

#Calculamos de nuevo el error absoluto medio
err_a = sum(abs(valist(1:length(dr_a), :)-dr_a))/length(dr_a);
err_b = sum(abs(vblist(1:length(dr_a), :)-dr_b))/length(dr_b);
err_c = sum(abs(vclist(1:length(dr_a), :)-dr_c))/length(dr_c);

fprintf("Error Absoluto Medio Derivada tras Usar Formula Extremo\n")
fprintf("Error A %f | %f | %f \n", err_a);
fprintf("Error B %f | %f | %f \n", err_b);
fprintf("Error C %f | %f | %f \n", err_c); #Lo hemos logrado reducir!!

#Usamos ahora la funcion cumtrapz para calcular la integral de v usando el metodo del trapezoide acumulativo
intv_a = cumtrapz(0:dt:(t-dt), valist);
intv_b = cumtrapz(0:dt:(t-dt), vblist);
intv_c = cumtrapz(0:dt:(t-dt), vclist);

#Calculamos el error absoluto medio
err_va = sum(abs(alist(1:length(intv_a)-1, :)-intv_a(2:length(intv_a), :)))/length(intv_a);
err_vb = sum(abs(blist(1:length(intv_a)-1, :)-intv_b(2:length(intv_a), :)))/length(intv_b);
err_vc = sum(abs(clist(1:length(intv_a)-1, :)-intv_c(2:length(intv_a), :)))/length(intv_c);

#Errores enormes y lo peor es que hay ejes donde los datos ni se parecen

fprintf("\nError Absoluto Medio Integral Metodo Trapezoide \n")
fprintf("Error A %f | %f | %f \n", err_va);
fprintf("Error B %f | %f | %f \n", err_vb);
fprintf("Error C %f | %f | %f \n", err_vc);

#Metodo extra Runge-Kutta

#Definimos nuestro sistema de ecuaciones diferenciales
function derivs = Ec_tres_cuerpos(t, w, G, m1, m2 ,m3)
  r1 = w(1:3);
  r2 = w(4:6);
  r3 = w(7:9);
  v1 = w(10:12);
  v2 = w(13:15);
  v3 = w(16:18);

  r12 = norm(r2-r1);
  r13 = norm(r3-r1);
  r23 = norm(r3-r2);

  dv1 = G*((m2*(r2-r1))/r12^3 + (m3*(r3-r1))/r13^3); #a = F/m
  dv2 = G*((m1*(r1-r2))/r12^3 + (m3*(r3-r2))/r23^3);
  dv3 = G*((m1*(r1-r3))/r13^3 + (m2*(r2-r3))/r23^3);

  dr1 = v1;
  dr2 = v2;
  dr3 = v3;

  derivs = [dr1; dr2; dr3; dv1; dv2; dv3];

endfunction

init_params = [r0_a, r0_b, r0_c, v0_a, v0_b, v0_c]; #Definimos los parametros iniciales
t_span = 0:dt:(t-dt);

#Usamos ode45 para que consiga las soluciones al sistema usando RK4
[~, sol] = ode45(@(t,w) Ec_tres_cuerpos(t, w, G, Ma, Mb, Mc), t_span, init_params);

r1_sol = sol(:, 1:3);
r2_sol = sol(:, 4:6);
r3_sol = sol(:, 7:9);

figure(5)
plot3(r1_sol(:, 1), r1_sol(:, 2), r1_sol(:, 3), 'b', r2_sol(:, 1), r2_sol(:, 2), r2_sol(:, 3), 'r', r3_sol(:, 1), r3_sol(:, 2), r3_sol(:, 3), 'g')
xlim([-1e12, 1e12])
ylim([-1.5e12, 1.5e12])
title("Trayectoria 3 cuerpos usando RK4");

figure(6)
hold on;
Animate(r1_sol, r2_sol, r3_sol);
title("Simulacion 3 Cuerpos RK4");

