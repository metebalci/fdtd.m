function fdtd_1d_run(layers, length_unit, d_min, f_max, update_freq)

    % universtal constants
    c0 = 299793458; % m/s
    
    % simulation constants
    % these can be modified if needed
    % or taken as input parameters
    N_wavelength = 20;
    N_dimension = 4;
    n_bc = 1;
    source_z = 2;
    spacer_region = 10;
    
    % delete any waitbar hanged from any previous session
    delete(findall(0,'type','figure','tag','TMWWaitbar'));
    wb = waitbar(0, "Initializing...", "Name", "Simulation Progress", "CreateCancelBtn", "setappdata(gcbf, 'cancel', 1)");
    
    % cache some quick calculations
    er_max = max(layers(:, 2));
    ur_max = max(layers(:, 3));
    n_max = sqrt(er_max * ur_max);

    % <--- COMPUTE GRID RESOLUTION --->
    lambda_min = c0 / n_max / f_max;

    d_c = d_min;

    dz_lambda = lambda_min / N_wavelength;
    dz_dmin = d_min / N_dimension;
    dz = min(dz_lambda, dz_dmin);
    % change dz to snap to critical dimension
    dz = d_c / ceil(d_c / dz);
    % check if d_c is integer multiple of dz
    assert (rem(d_c, dz) == 0);
    % <----->

    % <--- BUILD ENVIRONMENT --->
    % reflectance record, source inject, pre-space, post-space, tranmittance record
    Nz = 1 + 1 + spacer_region + spacer_region + 1;
    for m = 1:size(layers, 1)
        thickness = layers(m, 1)*length_unit;
        % check if start and finish are integer multiples of dz
        assert (rem(thickness, dz) == 0);
        Nz = Nz + thickness./dz;
    end
    assert (floor(Nz) == Nz, sprintf("Nz is not integer: %d", Nz));
    Nz = floor(Nz);
    ERyy = ones(1, Nz);
    URxx = ones(1, Nz);
    % epsilon relative
    pos = 3+spacer_region;
    for m = 1:size(layers, 1)
        thickness = floor(layers(m, 1)*length_unit./dz);
        er = layers(m, 2);
        ur = layers(m, 3);
        ERyy(pos:pos+thickness) = er;
        URxx(pos:pos+thickness) = ur;
        pos = pos + thickness;
    end
    % <----->

    % <--- COMPUTE TIME RESOLUTION --->
    dt = (n_bc * dz) / (2 * c0);
    % <----->

    % <--- COMPUTE TOTAL SIMULATION TIME --->
    tprop = n_max * Nz * dz / c0;    
    tau = 1/(pi*f_max);
    T = 12 * tau + 5 * tprop;    
    STEPS = ceil(T / dt);
    assert(tau / dt >= 10);
    % <------>

    % <--- COMPUTE SOURCE --->
    t0 = 6 * tau;
    gE = (0:STEPS-1) .* dt;
    gE = exp(- ((gE - t0) ./ tau) .^ 2);
    A = - sqrt(ERyy(source_z) / URxx(source_z));
    n_src = sqrt(ERyy(source_z) * URxx(source_z));
    gH = (0:STEPS-1) .* dt;
    delta_t = (n_src * dz)/(2 * c0) + dt/2;
    gH = A .* exp(- ((gH - t0 + delta_t) ./ tau) .^ 2);
    % <------>

    % <--- STARTING SIMULATION --->

    % INITIALIZE CONSTANTS, LOOP-INVARIANTS and VARIABLES
    % update coefficients
    mEy = (c0*dt) ./ ERyy ./ dz;
    mHx = (c0*dt) ./ URxx ./ dz;
    % fields
    Ey = zeros(1, Nz);
    Hx = zeros(1, Nz);
    % initialize boundary terms
    H1 = 0; 
    %H2 = 0; % this is first time initialized in the loop before use
    E1 = 0;
    E2 = 0;
    % initialize fourier parameters
    fourier_upper_limit = min(floor(0.5/dt), f_max);
    fourier_resolution = floor(1/dt/STEPS);
    n_freq = max(floor(fourier_upper_limit / fourier_resolution), 200);
    freq = linspace(0, fourier_upper_limit, n_freq);
    K = exp(-1i*2*pi*dt.*freq);
    % Kc caches K^T so it is not calculated in every iteration
    Kc = K;
    fref = zeros(1, n_freq);
    ftra = zeros(1, n_freq);
    fsrc = zeros(1, n_freq);

    % setup figure
    wavesplot = figure(1);
    clf(wavesplot);
    box on;
    hold on;

    % draw materials
    pos =3+spacer_region;
    colors = [0.5 0.5 0.5; 0.7 0.7 0.7];
    for m = 1:size(layers, 1)
        thickness = floor(layers(m, 1)*length_unit./dz);
        rectangle("Position", [pos, -2, thickness, 4], "FaceColor", colors(mod(m, 2)+1, :));
        pos = pos + thickness;
    end

    eplot = plot(Ey, "-b");
    hplot = plot(Hx, "-r");
    title("fields 2D");
    xlabel("z");
    ylabel("field");
    axis([1 Nz -2 2]);
    legend([eplot hplot], ["E", "H"], "location", "southoutside");

  %figure 1;
  %box on;
  %hold on;
  %view(3);
  %eplot3d = plot3(zeros(1, Nz), Ey, [1:Nz]);
  %hplot3d = plot3(Hx, zeros(1, Nz), [1:Nz]);
  %title("fields 3D")
  %xlabel("H");
  %ylabel("E");
  %zlabel("z");
  %axis([-2 2 -2 2 1 Nz]);
  %legend("E", "H");

    % STARTING SIMULATION

    waitbar(0, wb, "Starting...");

    for T = 1:STEPS

        % save for perfect boundary calculation
        H2=H1; H1=Hx(1);

        % H update
        Ey_shifted = circshift(Ey, -1);
        Ey_shifted(Nz) = E2; 
        Hx = Hx + mHx .* (Ey_shifted - Ey);

        % TF/SF correction
        Hx(source_z - 1) = Hx(source_z - 1) - mHx(source_z - 1) * gE(T);

        % save for perfect boundary calculation
        E2=E1; E1=Ey(Nz);

        % E update
        Hx_shifted = circshift(Hx, 1);
        Hx_shifted(1) = H2;
        Ey = Ey + mEy .* (Hx - Hx_shifted);

        % TF/SF correction
        Ey(source_z) = Ey(source_z) - mEy(source_z) * gH(T);

        % fourier update
        % 1 is reflectance measuring point
        fref = fref + Kc .* Ey(1);
        % Nz is transmittance measuring point
        ftra = ftra + Kc .* Ey(Nz);
        % to normalize
        fsrc = fsrc + Kc .* gE(T);
        % next Kc, to eliminate ^ calculation everytime
        Kc = Kc .* K;

        if rem(T, update_freq) == 0
            if getappdata(wb, "cancel")
                break
            end
            % update figure
            set(eplot, "YData", Ey);
            set(hplot, "YData", Hx);
            %set(eplot3d, "YData", Ey);
            %set(hplot3d, "XData", Hx);
            drawnow();
            waitbar(T/STEPS, wb, "Running...");
        end

    end

    delete(wb);

    s11 = abs(fref ./ fsrc);
    s21 = abs(ftra ./ fsrc);

    splot = figure(2);
    clf(splot);
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
        
end  
