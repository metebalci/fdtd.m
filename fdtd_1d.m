close all
clear

function fdtd = fdtd_init (layers, length_unit)

  fdtd = struct();
  fdtd.c0 = 299793458; % m/s
  fdtd.layers = layers
  fdtd.length_unit = length_unit;
  fdtd.N_wavelength = 20;
  fdtd.N_dimension = 4;
  fdtd.n_bc = 1;
  fdtd.source_z = 2;
  fdtd.spacer_region = 10;

endfunction

% minimum feature size, gcd of layer positions/thickness
function ret = fdtd_priv_d_min (fdtd)
  % gcd does not work for single number
  if size(fdtd.layers, 1) == 1
    ret = fdtd.layers(1, 1) * fdtd.length_unit;
  else
    % assuming the minimum length given is more than a picometer
    v = num2cell(arrayfun(@floor, fdtd.layers(:, 1) .* fdtd.length_unit .* 1e12));
    ret = gcd(v{:}) ./ 1e12;
  endif
endfunction

% critical feature size, same as minimum feature_size
function ret = fdtd_priv_d_c (fdtd)
  ret = fdtd_priv_d_min(fdtd);
endfunction

function ret = fdtd_priv_er_max (fdtd)
  ret = max(fdtd.layers(:, 2));
endfunction

function ret = fdtd_priv_ur_max (fdtd)
  ret = max(fdtd.layers(:, 3));
endfunction

% maximum refractive index
function ret = fdtd_priv_n_max (fdtd)
  ret = sqrt(fdtd_priv_er_max(fdtd) * fdtd_priv_ur_max(fdtd));
endfunction

function dz = fdtd_priv_compute_grid_resolution(fdtd, f_max)

  lambda_min = fdtd.c0/fdtd_priv_n_max(fdtd)/f_max
  d_min = fdtd_priv_d_min(fdtd)
  d_c = fdtd_priv_d_c(fdtd)
  dz_lambda = lambda_min / fdtd.N_wavelength
  dz_dmin = d_min / fdtd.N_dimension
  dz = min(dz_lambda, dz_dmin)
  % change dz to snap to critical dimension
  dz = d_c / ceil(d_c / dz)
  % check if d_c is integer multiple of dz
  assert (rem(d_c, dz) == 0);

endfunction

function [Nz, ERyy, URxx] = fdtd_priv_build_structure(fdtd, dz)

  length_unit = fdtd.length_unit;
  layers = fdtd.layers;

  % reflectance record, source inject, pre-space, post-space, tranmittance record
  Nz = 1 + 1 + fdtd.spacer_region + fdtd.spacer_region + 1
  for m = 1:size(layers, 1);
    thickness = layers(m, 1)*length_unit;
    % check if start and finish are integer multiples of dz
    assert (rem(thickness, dz) == 0);
    Nz = Nz + thickness./dz;
  endfor
  assert (floor(Nz) == Nz);
  Nz = floor(Nz);
  ERyy = ones(1, Nz);
  URxx = ones(1, Nz);
  % epsilon relative
  pos = 3+fdtd.spacer_region;
  for m = 1:size(layers, 1)
    thickness = floor(layers(m, 1)*length_unit./dz);
    er = layers(m, 2);
    ur = layers(m, 3);
    ERyy(pos:pos+thickness) = er;
    URxx(pos:pos+thickness) = ur;
    pos = pos + thickness;
  endfor

endfunction

function dt = fdtd_priv_compute_time_step(fdtd, dz)

  dt = (fdtd.n_bc * dz) / (2 * fdtd.c0);
  
endfunction

function [T tau tprop] = fdtd_priv_compute_total_time(fdtd, Nz, dz, f_max)

  tprop = fdtd_priv_n_max(fdtd) * Nz * dz / fdtd.c0;
  printf("tprop: %d\n", tprop);
  tau = 1/(pi*f_max);
  T = 12 * tau + 5 * tprop;

endfunction

function [gE gH] = fdtd_priv_compute_source(fdtd, dz, dt, tau, STEPS, ERyy, URxx)

  t0 = 6 * tau
  gE = [0:STEPS-1] .* dt;
  gE = exp(- ((gE - t0) ./ tau) .^ 2);
  A = - sqrt(ERyy(fdtd.source_z) / URxx(fdtd.source_z));
  n_src = sqrt(ERyy(fdtd.source_z) * URxx(fdtd.source_z));
  gH = [0:STEPS-1] .* dt;
  delta_t = (n_src * dz)/(2 * fdtd.c0) + dt/2;
  printf("delta_t: %d\n", delta_t);
  gH = A .* exp(- ((gH - t0 + delta_t) ./ tau) .^ 2);

endfunction

function fdtd = fdtd_run(layers, length_unit, f_max)

  fdtd = fdtd_init(layers, length_unit)

  er_min = min(fdtd.layers(:, 2));
  er_max = max(fdtd.layers(:, 2));

  % compute grid resolution
  dz = fdtd_priv_compute_grid_resolution(fdtd, f_max);

  if (dz > 1e-3)
    printf("dz: %d mm\n", dz*1e3);
  elseif (dz > 1e-6)
    printf("dz: %d um\n", dz*1e6);
  else
    printf("dz: %d nm\n", dz*1e9);
  endif

  % build environment
  [Nz, ERyy, URxx] = fdtd_priv_build_structure(fdtd, dz);

  printf("Nz: %d\n", Nz);

  % compute time resolution
  dt = fdtd_priv_compute_time_step(fdtd, dz)

  if (dt > 1e-9)
    printf("dt: %d ns\n", dt*1e9);
  elseif (dt > 1e-12)
    printf("dt: %d ps\n", dt*1e12);
  else 
    printf("dt: %d fs\n", dt*1e15);
  endif

  % compute total simulation time
  [T tau tprop] = fdtd_priv_compute_total_time(fdtd, Nz, dz, f_max);
  printf("T: %d\n", T);
  printf("tau: %d\n", tau);
  STEPS = ceil(T / dt);
  printf("STEPS: %d\n", STEPS);

  assert(tau / dt >= 10);

  % compute source
  [gE gH] = fdtd_priv_compute_source(fdtd, dz, dt, tau, STEPS, ERyy, URxx);

  % STARTING SIMULATION

  % INITIALIZE CONSTANTS AND INITIAL E, H VALUES
  % update coefficients
  source_z = fdtd.source_z
  mEy = (fdtd.c0*dt) ./ ERyy;
  mHx = (fdtd.c0*dt) ./ URxx;
  mEy = mEy ./ dz;
  mHx = mHx ./ dz;
  % fields
  Ey = zeros(1, Nz);
  Hx = zeros(1, Nz);
  % initialize boundary terms
  H1 = H2 = 0;
  E1 = E2 = 0;
  % initialize fourier parameters
  fourier_upper_limit = min(floor(0.5/dt), f_max);
  printf("fourier_upper_limit: %d GHz\n", fourier_upper_limit/1e9);
  fourier_resolution = floor(1/dt/STEPS);
  printf("fourier_resolution: %d MHz\n", fourier_resolution/1e6);
  fourier_number_of_points = floor(fourier_upper_limit / fourier_resolution)
  printf("fourier_number_of_points: %d\n", fourier_number_of_points)
  n_freq = max(floor(fourier_upper_limit / fourier_resolution), 200);
  freq = linspace(0, fourier_upper_limit, n_freq);
  K = exp(-i*2*pi*dt.*freq);
  fref = zeros(1, n_freq);
  ftra = zeros(1, n_freq);
  fsrc = zeros(1, n_freq);

  % setup figure
  figure 1;
  box on;
  hold on;

  % draw materials
  pos =3+fdtd.spacer_region;
  for m = 1:size(fdtd.layers, 1)
    thickness = floor(fdtd.layers(m, 1)*fdtd.length_unit./dz);
    er = floor(fdtd.layers(m, 2));
    color = ones(1, 3) .* 0.4 + round((er-er_min)/(er_max-er_min)) .* 0.3;
    rectangle("Position", [pos, 2, thickness, -4], "FaceColor", color);
    pos = pos + thickness;
  endfor

  eplot = plot(Ey, "-b");
  hplot = plot(Hx, "-r");
  title("fields");
  xlabel("z");
  ylabel("field");
  axis([1 Nz -2 2]);
  legend([eplot hplot], {"E", "H"}, "location", "southoutside");

  % SIMULATION

  % Kc caches K^T so it is not calculated in every iteration
  Kc = K;

  update_freq = 100;
  t0 = clock();
  for T = 1:STEPS

    % save for perfect boundary calculation
    H2=H1; H1=Hx(1);

    %sequential implementation
    %for k=1:Nz-1
    %  Hx(k) = Hx(k) + mHx(k)*(Ey(k+1) - Ey(k));
    %endfor
    %perfect boundary
    %Hx(Nz) = Hx(Nz) + mHx(Nz)*(E2 - Ey(Nz));

    %vectoral implementation
    Ey_shifted = shift(Ey, -1);
    Ey_shifted(Nz) = E2; 
    Hx = Hx + mHx .* (Ey_shifted - Ey);

    % TF/SF correction
    Hx(source_z - 1) = Hx(source_z - 1) - mHx(source_z - 1) * gE(T);

    % save for perfect boundary calculation
    E2=E1; E1=Ey(Nz);

    %sequential implementation
    %perfect boundary
    %Ey(1) = Ey(1) + mEy(1)*(Hx(1) - H2);
    %for k=2:Nz
    %  Ey(k) = Ey(k) + mEy(k)*(Hx(k) - Hx(k-1));
    %endfor

    %vectoral implementation
    Hx_shifted = shift(Hx, 1);
    Hx_shifted(1) = H2;
    Ey = Ey + mEy .* (Hx - Hx_shifted);

    % TF/SF correction
    Ey(source_z) = Ey(source_z) - mEy(source_z) * gH(T);

    % update fourier transforms
    % 1 is reflectance measuring point
    fref = fref + Kc .* Ey(1);
    % Nz is transmittance measuring point
    ftra = ftra + Kc .* Ey(Nz);
    % to normalize
    fsrc = fsrc + Kc .* gE(T);
    % next Kc, to eliminate ^ calculation everytime
    Kc = Kc .* K;

    if rem(T, update_freq) == 0
      % update figure
      set(eplot, "YData", Ey);
      set(hplot, "YData", Hx);
      drawnow();
      throughput = floor(1/(etime(clock(), t0)/update_freq));
      remaining_seconds = floor((STEPS-T)/throughput);
      printf("T: %d (%d steps/sec, %d seconds remained)\n", T, throughput, remaining_seconds)
      t0 = clock();
    end

  endfor

  s11 = abs(fref ./ fsrc);
  s21 = abs(ftra ./ fsrc);

  figure 2;
  box on;
  hold on;
  plot(freq, 20*log(s11), "-r");
  plot(freq, 20*log(s21), "-b");
  plot(freq, s11+s21, "-k");
  title("s-parameters");
  xlabel("f");
  ylabel("s (dB)");
  axis([0 fourier_upper_limit -100 20]);
  legend("s11", "s21", "conservation (s11+s21)", "location", "southoutside");
  
endfunction

% example here: https://empossible.net/wp-content/uploads/2020/01/Lecture-Review-Walkthrough-of-1D-FDTD.pdf
function fdtd_demo1()

  close all;

  fdtd_run([
    304800 2 6
    ], 1e-6, 1e9)

endfunction

% anti-reflection example #1 here: https://empossible.net/wp-content/uploads/2020/01/Lecture-Examples-of-1D-FDTD.pdf
function fdtd_demo2()

  close all;

  c0 = 299793458; % m/s

  e_r = 12
  e_air = 1
  e_ar = sqrt(e_r * e_air)
  n_ar = sqrt(e_ar)
  f0 = 2.4e9;
  lambda_0 = c0 / f0
  d_ar = floor(lambda_0 / (4 * n_ar) / 1e-3)

  fdtd_run([
    d_ar e_ar 1
    304.8 e_r 1
    d_ar e_ar 1
    ], 1e-3, 4e9)

endfunction

% anti-reflection example #2 here: https://empossible.net/wp-content/uploads/2020/01/Lecture-Examples-of-1D-FDTD.pdf
function fdtd_demo3()

  close all;

  c0 = 299793458; % m/s

  e_sin = 1.5*1.5;
  e_sio2 = 2*2;
  d1 = 163
  d2 = 122

  fdtd_run([
    d1 e_sin 1
    d2 e_sio2 1
    d1 e_sin 1
    d2 e_sio2 1
    d1 e_sin 1
    d2 e_sio2 1
    d1 e_sin 1
    d2 e_sio2 1
    d1 e_sin 1
    d2 e_sio2 1
    d1 e_sin 1
    d2 e_sio2 1
    d1 e_sin 1
    d2 e_sio2 1
    d1 e_sin 1
    d2 e_sio2 1
    d1 e_sin 1
    d2 e_sio2 1
    d1 e_sin 1
    d2 e_sio2 1
    ], 1e-9, 2*c0/980e-9)

endfunction
