
function E = computeHF_E(Kc, L, Rays, v)

    % Function to compute E based on the given equation
    % Inputs:
    % Kc   - Stiffness coefficient
    % L    - Characteristic length
    % Rays - Number of rays
    % v    - Poisson's ratio
    % Output:
    % E    - Computed value

    % Compute E
    E = ((Kc / L) * (Rays / (2 * pi))) * (1 + v) * (1 - 2 * v);
end