clear; clc; close all;

%% Settings
load('motor_flux_inductance_maps_sampled.mat');   % map_*_sample, params
compare_with_real = true;     % also compute MTPA from the full-resolution map
real_file         = 'motor_flux_inductance_maps_real.mat';
save_lut          = true;     % save torque -> (id, iq, psi_s) table
N_is   = 100;                 % number of current amplitudes
N_beta = 7201;                % angle points per amplitude (0.05 deg)
export = 1;
Udc    = 800;                 % DC link voltage [V]
Us_max = Udc/sqrt(3);         % max phase voltage amplitude [V]
Imax   = 15;                  % current limit (amplitude) [A]
Imap   = 16;                  % maps cover +-Imap A
p      = params.p;
plttrq = 1;
%% Unpack
id_s   = map_PSID_sample(1, :, 1);
iq_s   = map_PSID_sample(:, 1, 2);
PSID_s = map_PSID_sample(:, :, 3);
PSIQ_s = map_PSIQ_sample(:, :, 3);

%% MTPA from sampled maps
R = mtpa_search(id_s, iq_s, PSID_s, PSIQ_s, p, Imax, N_is, N_beta, plttrq, export);

%% Torque-to-current look-up table along the MTPA
T_grid = linspace(0, R.T(end), 64)';
[T_u, iu] = unique(R.T);
if numel(T_u) < numel(R.T)
    warning('MTPA torque not strictly increasing - LUT built from unique points.');
end
LUT.T    = T_grid;
LUT.id   = interp1(T_u, R.id(iu),   T_grid, 'pchip');
LUT.iq   = interp1(T_u, R.iq(iu),   T_grid, 'pchip');
LUT.psis = interp1(T_u, R.psis(iu), T_grid, 'pchip');
if save_lut
    save('mtpa_lut.mat', 'LUT');
end

%% Fine grid of torque for plotting
id_f = linspace(-Imap, Imap, 321);
iq_f = linspace(-Imap, Imap, 321);
[IDf, IQf] = meshgrid(id_f, iq_f);
PDf = interp2(id_s, iq_s, PSID_s, IDf, IQf, 'spline');
PQf = interp2(id_s, iq_s, PSIQ_s, IDf, IQf, 'spline');
Tf  = 1.5*p*(PDf.*IQf - PQf.*IDf);

%% Plot dq currents with clean torque contours and accurate legend
figure('Name','MTPA and Torque Contours','Color','w'); 
hold on; grid on; box on;

[C, h_c] = contour(IDf, IQf, Tf, 15, 'Color', [0.6 0.6 0.6], 'LineWidth', 1);

h_dummy = plot(nan, nan, 'Color', [0.6 0.6 0.6], 'LineWidth', 1);

thc = linspace(0, 2*pi, 400);

h_imax = plot(Imax*cos(thc), Imax*sin(thc), 'Color', [0 0.45 0.74], 'LineWidth', 2);

h_mtpa = plot(R.id, R.iq, 'Color', [0.85 0.33 0.1], 'LineWidth', 2.5);

axis equal; 
xlim([-Imap Imap]); 
ylim([-Imap Imap]);
xlabel('i_d [A]'); 
ylabel('i_q [A]');
legend([h_dummy, h_imax, h_mtpa], ...
       {'Torque Contours', sprintf('I_{max} Limit (%d A)', Imax), 'MTPA Trajectory'}, ...
       'Location', 'northeastoutside');
if export == 1
    matlab2tikz('dq_currents_plot.tex', 'standalone', true);
end

%% Local function: brute-force MTPA search on interpolated maps
function R = mtpa_search(id_g, iq_g, PSID, PSIQ, p, Imax, N_is, N_beta, plttrq, export)
    is_vec = linspace(0, Imax, N_is)';
    beta   = linspace(0, 2*pi, N_beta);
    ID_q = is_vec*cos(beta);
    IQ_q = is_vec*sin(beta);
    PD = interp2(id_g, iq_g, PSID, ID_q, IQ_q, 'spline');
    PQ = interp2(id_g, iq_g, PSIQ, ID_q, IQ_q, 'spline');
    Tq = 1.5*p*(PD.*IQ_q - PQ.*ID_q);
    [Tm, kb] = max(Tq, [], 2);
    lin = sub2ind(size(Tq), (1:N_is)', kb);
    R.is   = is_vec;
    R.T    = Tm;
    R.id   = ID_q(lin);
    R.iq   = IQ_q(lin);
    R.psid = PD(lin);
    R.psiq = PQ(lin);
    R.psis = hypot(R.psid, R.psiq);
    R.beta = atan2d(R.iq, R.id);
    R.kt   = Tm ./ is_vec;
    if plttrq == 1
        figure('Name','Electromagnetic Torque Map','Color','w');
        surf(ID_q, IQ_q, Tq, 'EdgeColor','none'); 
        colormap(hot);
        colorbar;
        hold on;
        
        % 2. Thin black dashed line on top for contrast against light areas
        plot3(R.id, R.iq, R.T, 'Color', [0.8 0 0.8], 'LineWidth', 3);
        xlabel('i_d [A]'); 
        ylabel('i_q [A]'); 
        zlabel('Torque [Nm]');
        view(40,30); 
        grid on;
        if export == 1
            matlab2tikz('Torque_map.tex', 'standalone', true);
        end
    end
end