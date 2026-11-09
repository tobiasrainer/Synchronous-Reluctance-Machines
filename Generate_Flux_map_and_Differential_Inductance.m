clear; clc; close all;

plt = 0;
export = 0;
save_data = 0;

%% Init parameters
params.p      = 2;       % pole pairs
params.seed   = 1;       % random seed
params.Ldinf  = 0.010;   % d inductance far in saturation [H]
params.Lqinf  = 0.018;   % q inductance far in saturation [H]
params.Lqpeak = 0.047;   % target peak of Lqq [H]
params.wd     = 0.7;     % Fluxmap parameter
params.K      = 8;       % Fluxmap parameter

Udc    = 800;                 % DC link voltage [V]
Us_max = Udc/sqrt(3);         % max phase voltage amplitude [V]
Imax   = 15;                  % current limit (amplitude) [A]
Imap   = 16;                  % maps cover +-Imap A
n_rpm  = 9000;                % nominal speeds [rpm]
n_list = [4000 6000 9000 12000];     % speeds for the voltage-limit curves [rpm]
n_env  = linspace(200, 15000, 150);  % speeds for torque-speed envelope [rpm]

%% Random parameters for flux model   
rng(params.seed);
K  = params.K;  Kh = K/2;
th = [pi/2 + (2*rand(1,Kh)-1)*0.5, (2*rand(1,Kh)-1)*0.5];  % ridge directions
m.a    = cos(th);
m.b    = sin(th);
m.s    = 6 + 8*rand(1,K);                  % ridge widths [A]
m.ci   = (2*rand(1,K)-1)*8;                % ridge offsets in id [A]
m.cq   = (2*rand(1,K)-1)*8;                % ridge offsets in iq [A]
m.w    = (0.5 + rand(1,K)) .* [ones(1,Kh), params.wd*ones(1,Kh)];
m.Ldinf = params.Ldinf;
m.Lqinf = params.Lqinf;
m.g    = 1;

N    = 321;                            % odd -> id = 0 / iq = 0 lies on the grid
id_v = linspace(-Imap, Imap, N);
iq_v = linspace(-Imap, Imap, N);
[ID, IQ] = meshgrid(id_v, iq_v);       % rows = iq, columns = id
Icur = (ID.^2 + IQ.^2) <= Imax^2;      % inside current circle

[~, ~, ~, Lqq1] = fluxmodel(ID, IQ, m);                % raw maps with g = 1
m.g = (params.Lqpeak - m.Lqinf) / (max(Lqq1(Icur)) - m.Lqinf);
[PSID, PSIQ, LDD, LQQ, LDQ] = fluxmodel(ID, IQ, m);    % final maps

PSIS2 = PSID.^2 + PSIQ.^2;
T     = 1.5*params.p*(PSID.*IQ - PSIQ.*ID);               % torque map [Nm]

ix = 1:20:N;

%% Plots
if plt == 1
    figure('Name','Differential Inductance L_{dd}','Color','w');
    surf(ID(ix,ix), IQ(ix,ix), LDD(ix,ix)); 
    colormap(hot); 
    colorbar;
    xlabel('i_d [A]'); 
    ylabel('i_q [A]'); 
    zlabel('L_{dd} [H]');
    view(40,30); 
    grid on;
    if export == 1
        matlab2tikz('Ldd_plot.tex', 'standalone', true);
    end

    figure('Name','Differential Inductance L_{qq}','Color','w');
    surf(ID(ix,ix), IQ(ix,ix), LQQ(ix,ix)); 
    colormap(hot); 
    colorbar;
    xlabel('i_d [A]'); 
    ylabel('i_q [A]'); 
    zlabel('L_{qq} [H]');
    view(40,30); 
    grid on;
    if export == 1
        matlab2tikz('Lqq_plot.tex', 'standalone', true);
    end

    figure('Name','Differential Inductance L_{dq}','Color','w');
    surf(ID(ix,ix), IQ(ix,ix), LDQ(ix,ix)); 
    colormap(hot); 
    colorbar;
    xlabel('i_d [A]'); 
    ylabel('i_q [A]'); 
    zlabel('L_{dq} = L_{qd} [H]');
    view(40,30); 
    grid on;
    if export == 1
        matlab2tikz('Ldq_plot.tex', 'standalone', true);
    end

    figure('Name','Flux Linkage \psi_d','Color','w');
    surf(ID, IQ, PSID, 'EdgeColor','none'); 
    colormap(hot); 
    colorbar;
    xlabel('i_d [A]'); 
    ylabel('i_q [A]'); 
    zlabel('\psi_d [Vs]');
    view(40,30); 
    grid on;
    if export == 1
        matlab2tikz('PSID_plot.tex', 'standalone', true);
    end

    figure('Name','Flux Linkage \psi_q','Color','w');
    surf(ID, IQ, PSIQ, 'EdgeColor','none'); 
    colormap(hot); 
    colorbar;
    xlabel('i_d [A]'); 
    ylabel('i_q [A]'); 
    zlabel('\psi_q [Vs]');
    view(40,30); 
    grid on;
    if export == 1
        matlab2tikz('PSIQ_plot.tex', 'standalone', true);
    end
end

%% Save all map matrices to .mat file
if save_data == 1
    save_motor_maps('motor_flux_inductance_maps_real.mat','motor_flux_inductance_maps_sampled.mat', ID, IQ, LDD, LQQ, LDQ, PSID, PSIQ, ix, params);
end

%% Local function: flux linkages and differential inductances
function [psid, psiq, Ldd, Lqq, Ldq] = fluxmodel(id, iq, m)
    psid = m.Ldinf*id;
    psiq = m.Lqinf*iq;
    Ldd  = m.Ldinf*ones(size(id));
    Lqq  = m.Lqinf*ones(size(id));
    Ldq  = zeros(size(id));
    for k = 1:numel(m.a)
        u  = (m.a(k)*(id - m.ci(k)) + m.b(k)*(iq - m.cq(k)))/m.s(k);
        u0 = (-m.a(k)*m.ci(k) - m.b(k)*m.cq(k))/m.s(k);              % u at (0,0)
        f  = m.g*m.w(k)*m.s(k)*sqrt(pi)/2 .* (erf(u) - erf(u0));     % shifted so psi(0,0) = 0
        psid = psid + f*m.a(k);
        psiq = psiq + f*m.b(k);
        e  = m.g*m.w(k)*exp(-u.^2);
        Ldd = Ldd + e*m.a(k)^2;
        Lqq = Lqq + e*m.b(k)^2;
        Ldq = Ldq + e*m.a(k)*m.b(k);
    end
end

%% Local function: generate and save all map structures
function save_motor_maps(filename_1, filename_2, ID, IQ, LDD, LQQ, LDQ, PSID, PSIQ, ix, params)
    map_LDD_full = cat(3, ID, IQ, LDD);
    map_LQQ_full = cat(3, ID, IQ, LQQ);
    map_LDQ_full = cat(3, ID, IQ, LDQ);
    map_PSID_full = cat(3, ID, IQ, PSID);
    map_PSIQ_full = cat(3, ID, IQ, PSIQ);
    
    ID_sample   = ID(ix, ix);
    IQ_sample   = IQ(ix, ix);
    
    map_LDD_sample  = cat(3, ID_sample, IQ_sample, LDD(ix, ix));
    map_LQQ_sample  = cat(3, ID_sample, IQ_sample, LQQ(ix, ix));
    map_LDQ_sample  = cat(3, ID_sample, IQ_sample, LDQ(ix, ix));
    map_PSID_sample = cat(3, ID_sample, IQ_sample, PSID(ix, ix));
    map_PSIQ_sample = cat(3, ID_sample, IQ_sample, PSIQ(ix, ix));
    
    save(filename_1, ...
         'map_LDD_full', 'map_LQQ_full', 'map_LDQ_full', ...
         'map_PSID_full', 'map_PSIQ_full', ...
         'params', '-v7.3');
    save(filename_2, ...
         'map_LDD_sample', 'map_LQQ_sample', 'map_LDQ_sample', ...
         'map_PSID_sample', 'map_PSIQ_sample', ...
         'params', '-v7.3');
end