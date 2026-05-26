function custom_light_cone_projection
    % custom_light_cone_projection.m
    % Calculates the projection of a light cone onto a slanted screen.
    
    clc; clear; close all;

    % ==========================================
    % 1. Custom Input Parameters
    % ==========================================
    L0_target = 10;         % Point-to-point distance L0 (meters)
    theta0_deg = 45;        % Grazing angle of the central axis (degrees)
    Omega = 0.63;           % Solid angle of the light cone (steradians)
    
    % ==========================================
    % 2. Geometry Conversions
    % ==========================================
    phi = acos(1 - Omega / (2*pi)); % Calculate cone half-angle
    theta0 = deg2rad(theta0_deg);
    
    fprintf('=== Input Parameters ===\n');
    fprintf('Distance to center L0 = %.2f m\n', L0_target);
    fprintf('Grazing angle theta_0 = %.2f deg\n', theta0_deg);
    fprintf('Cone half-angle phi = %.2f deg\n\n', rad2deg(phi));
    
    if theta0 <= phi
        error('Grazing angle is smaller than the cone half-angle. The projection is not a closed ellipse!');
    end

    % ==========================================
    % 3. Core Calculations
    % ==========================================
    alpha_deg = linspace(0, 360, 361); % Azimuthal angle on the screen
    alpha = deg2rad(alpha_deg);
    
    % Solve the quadratic equation A*R^2 + B*R + C = 0 for the edge distance R
    A = (cos(alpha).^2) .* (cos(theta0)^2) - cos(phi)^2;
    B = 2 * L0_target .* cos(alpha) .* cos(theta0) .* (sin(phi)^2);
    C = L0_target^2 * sin(phi)^2;
    
    % Extract the positive root for the edge distance
    R = (-B - sqrt(B.^2 - 4.*A.*C)) ./ (2.*A);
    
    % Calculate XY physical projection coordinates (Origin is at the optical center O)
    X = R .* cos(alpha);
    Y = R .* sin(alpha);
    
    % Calculate the incidence angle (psi) at each edge point
    H_source = L0_target * sin(theta0); % Vertical height of the source relative to the screen
    % Slant distance D from source to edge point
    D = sqrt(R.^2 + L0_target^2 + 2 .* R .* L0_target .* cos(alpha) .* cos(theta0));
    psi_deg = rad2deg(asin(H_source ./ D));

    % ==========================================
    % 4. Console Outputs (Macroscopic Parameters)
    % ==========================================
    idx_far = 1;     % 0 degrees: Far end
    idx_near = 181;  % 180 degrees: Near end
    idx_side = 91;   % 90 degrees: Side
    
    fprintf('=== Point-Specific Results ===\n');
    fprintf('[Far End]  Distance R: %.3f m, Incidence Angle: %.2f deg\n', R(idx_far), psi_deg(idx_far));
    fprintf('[Near End] Distance R: %.3f m, Incidence Angle: %.2f deg\n', R(idx_near), psi_deg(idx_near));
    fprintf('[Side Edge] Distance R: %.3f m, Incidence Angle: %.2f deg\n', R(idx_side), psi_deg(idx_side));
    
    % Calculate true geometric center and axes of the ellipse
    center_x = (X(idx_far) + X(idx_near)) / 2;
    length_major = X(idx_far) - X(idx_near);
    length_minor = max(Y) * 2;
    
    fprintf('\n=== Macroscopic Footprint ===\n');
    fprintf('Total Major Axis Length (X-direction): %.3f m\n', length_major);
    fprintf('Total Minor Axis Length (Y-direction): %.3f m\n', length_minor);
    fprintf('Geometric Center Offset (X-shift): %.3f m\n', center_x);

    % ==========================================
    % 5. Visualization
    % ==========================================
    figure('Color', 'w', 'Position', [100, 100, 1400, 450]);
    
    % --- Subplot 1: XY Physical Footprint ---
    subplot(1, 3, 1);
    plot(X, Y, 'b-', 'LineWidth', 2); hold on;
    plot(0, 0, 'k+', 'MarkerSize', 10, 'LineWidth', 2); % Optical Center O
    plot(center_x, 0, 'r.', 'MarkerSize', 15);          % Geometric Center
    
    text(X(idx_far), 0, ' Far Edge', 'VerticalAlignment', 'middle');
    text(X(idx_near), 0, 'Near Edge ', 'HorizontalAlignment', 'right');
    
    grid on; axis equal;
    xlabel('X Axis (m) [Slant Direction]'); 
    ylabel('Y Axis (m) [Transverse]');
    title('XY Footprint on the Screen');
    legend('Boundary', 'Optical Center (O)', 'Geometric Center', 'Location', 'best');
    
    % --- Subplot 2: Edge Distance R vs Alpha ---
    subplot(1, 3, 2);
    plot(alpha_deg, R, 'g-', 'LineWidth', 2);
    
    grid on;
    xlabel('Azimuthal Angle \alpha (degrees)');
    ylabel('Distance from Center R (m)');
    title('Edge Distance R vs. Angle \alpha');
    xlim([0 360]);
    xticks(0:90:360);
    xticklabels({'0 (Far)', '90 (Side)', '180 (Near)', '270 (Side)', '360'});
    
    % --- Subplot 3: Incidence Angle Psi vs Alpha ---
    subplot(1, 3, 3);
    plot(alpha_deg, psi_deg, 'r-', 'LineWidth', 2); hold on;
    yline(theta0_deg, 'k:', 'Reference Grazing Angle \theta_0');
    
    grid on;
    xlabel('Azimuthal Angle \alpha (degrees)');
    ylabel('Incidence Angle \psi (degrees)');
    title('Incidence Angle vs. Angle \alpha');
    xlim([0 360]);
    xticks(0:90:360);
    xticklabels({'0 (Far)', '90 (Side)', '180 (Near)', '270 (Side)', '360'});
end