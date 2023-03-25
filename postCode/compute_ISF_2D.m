function ISF = compute_ISF_2D(x1, y1, x2, y2, k, Ly)
% This function computes the intermediate scattering function for a 2D
% system, given the positions of particles x and y, using wavevectors in 8
% directions.

% Define parameters
N = length(x1);  % Number of particles

% Compute wavevectors in 8 directions
k_r=k;
k_t = [0, pi/4, pi/2, 3*pi/4, pi, 5*pi/4, 3*pi/2, 7*pi/4];
kx = k_r*cos(k_t);
ky = k_r*sin(k_t);

%k = [kx, ky];

% Compute distances between all pairs of particles
dx = x1 - x2;
dy = y1 - y2;
dy = mod(dy + Ly/2, Ly) - Ly/2;
r = sqrt(dx.^2 + dy.^2);

% Compute the structure factor for each wavevector
S = zeros(8,1);
for ii = 1:8
    S(ii) = sum(exp(1i*(kx(ii)*dx + ky(ii)*dy)),'all') / N;
end

% Compute the intermediate scattering function
ISF = sum(abs(S).^2)/8;

end
