%constantes
w_0 = 2;
zeta = 1;
F = 4;
w = 3;
m = 1;
%limites e size-step
a = 0;
b = 20;
h = 0.1;
n = (b-a)/h;
%cond. iniciais
x0 = 5;
v0 = 0;
t=[a zeros(1,n)];
x=[x0 zeros(1,n)];
v=[v0 zeros(1,n)];
%funções para o método
f1 = @(t,x,v) (v); %dx/dt = v
f2 = @(t,x,v) ((1/m)*F*sin(w*t)-2*w_0*zeta*v-(w_0^2)*x);

%-----------------------método--------------------------
for i = 1:n+1
    t(i+1) = t(i) + h;
    %determinaçãos dos k's
    k1_x = f1(t(i),x(i),v(i));
    k1_v = f2(t(i),x(i),v(i));
    k2_x = f1(t(i)+h,x(i)+k1_x*h,v(i)+k1_v*h);
    k2_v = f2(t(i)+h,x(i)+k1_x*h,v(i)+k1_v*h);
    %determinação dos x's e v's seguintes
    x(i+1) = x(i) + (0.5*k1_x+0.5*k2_x)*h;
    v(i+1) = v(i) + (0.5*k1_v+0.5*k2_v)*h;
    %imprimir valores
    fprintf('%5.4f  %11.8f %11.8f %11.8f %11.8f %11.8f %11.8f\n',...
        t(i), x(i), v(i), k1_x, k1_v, k2_x, k2_v);
    
end


%desenhar gráfico x(t)

plot(t,x); grid on;
title("Oscilador forçado com atrito, com m = 1kg");
hold on;
%desenhar gráfico v(t)

plot(t,v); grid on;
hold off;
legend('x','v');