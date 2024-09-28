function result = cel(kc, p, c, s)
    if kc == 0
        result = NaN;
        return;
    end

    errtol = 1e-6;
    k = abs(kc);
    pp = p;
    cc = c;
    ss = s;
    em = 1;

    if p > 0
        pp = sqrt(p);
        ss = s / pp;
    else
        f = kc^2;
        q = 1 - f;
        g = 1 - pp;
        f = f - pp;
        q = q * (ss - c * pp);
        pp = sqrt(f / g);
        cc = (c - ss) / g;
        ss = -q / (g^2 * pp) + cc * pp;
    end

    f = cc;
    cc = cc + ss / pp;
    g = k / pp;
    ss = 2 * (ss + f * g);
    pp = g + pp;
    g = em;
    em = k + em;
    kk = k;

    while abs(g - k) > g * errtol
        k = 2 * sqrt(kk);
        kk = k * em;
        f = cc;
        cc = cc + ss / pp;
        g = kk / pp;
        ss = 2 * (ss + f * g);
        pp = g + pp;
        g = em;
        em = k + em;
    end

    result = (pi / 2) * -(ss + cc * em) / (em * (em + pp));
end

% Parámetros constantes
mu0 = 4 * pi * 1e-7;  % T·m/A (permeabilidad del vacío)
n = 28;  % espiras/m
I = 0.7;  % corriente en A (2bnl)
b = 6e-3;  % longitud del solenoide (m)
rho = 5*1e-3;  % posición radial (m)
Z = 2.15*1e-2;  % posición axial (m)

% Calcular B_0
B0 = (mu0 / pi) * n * I;

% Rango de a: a = 3.5*10^{-3} + 3*10^{-4}*a (a desde 0 hasta 20)
a_vals = 3.5e-3 + (0:18) * 3e-4;

% Inicializar matrices para almacenar resultados de B_rho, B_z y B_total
B_rho_vals = zeros(1, length(a_vals));
B_z_vals = zeros(1, length(a_vals));
B_total_vals = zeros(1, length(a_vals));

for i = 1:length(a_vals)
    a = a_vals(i);
    
    % Calcular Z_+ y Z_-
    Z_plus = Z + b;
    Z_minus = Z - b;
    
    % Calcular k_+ y k_-
    k_plus = sqrt((Z_plus^2 + (a - rho)^2) / (Z_plus^2 + (a + rho)^2));
    k_minus = sqrt((Z_minus^2 + (a - rho)^2) / (Z_minus^2 + (a + rho)^2));
    
    % Calcular α_+ y α_-
    alpha_plus = a / sqrt(Z_plus^2 + (rho + a)^2);
    alpha_minus = a / sqrt(Z_minus^2 + (rho + a)^2);
    
    % Calcular β_+ y β_-
    beta_plus = Z_plus / (Z_plus^2 + (rho + a)^2);
    beta_minus = Z_minus / (Z_minus^2 + (rho + a)^2);
    
    % Calcular γ
    gamma = (a - rho) / (a + rho);
    
    % Calcular la componente radial B_rho
    B_rho_vals(i) = B0 * (alpha_plus * cel(k_plus, 1, 1, -1) - alpha_minus * cel(k_minus, 1, 1, -1));
    
    % Calcular la componente longitudinal B_z
    B_z_vals(i) = (B0 * a / (a + rho)) * (beta_plus * cel(k_plus, gamma^2, 1, gamma) - beta_minus * cel(k_minus, gamma^2, 1, gamma));
    
    % Calcular la magnitud total del campo magnético B
    B_total_vals(i) = sqrt(B_rho_vals(i)^2 + B_z_vals(i)^2);
end

% Calcular las sumas de las componentes B_rho, B_z y B_total
sum_B_rho = sum(B_rho_vals);
sum_B_z = sum(B_z_vals);
sum_B_total = sum(B_total_vals);

% Mostrar resultados de las sumas
fprintf('Suma de B_rho: %.6e T\n', sum_B_rho);
fprintf('Suma de B_z: %.6e T\n', sum_B_z);
fprintf('Suma de B_total: %.6e T\n', sum_B_total);

% Mostrar resultados de los valores individuales
fprintf('Resultados para B_rho, B_z y B_total:\n');
for i = 1:length(a_vals)
    fprintf('a = %.4e m: B_rho = %.6e T, B_z = %.6e T, B_total = %.6e T\n', a_vals(i), B_rho_vals(i), B_z_vals(i), B_total_vals(i));
end


subplot(3, 1, 3);
plot(a_vals, B_total_vals, '-o');
xlabel('a (m)');
ylabel('B_{total} (T)');
title('Magnitud total del campo magnético B_{total} en función de a');

