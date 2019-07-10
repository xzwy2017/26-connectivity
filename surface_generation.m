function [z , PixelWidth , PSD] = surface_generation(sigma, H, Lx, m , n,qs, qr,qL)
% Generates artifial randomly rough surfaces with given parameters.
% The code is adapted from https://ch.mathworks.com/matlabcentral/fileexchange/60817-surface-generator-artificial-randomly-rough-surfaces
% Copyright (c) 2016, Mona Mahboob Kanafi
% All rights reserved.

% ======================= inputs
% parameters (in SI units)
% sigma: standard deviation , i.e. root-mean-square roughness Rq(m)
% H: Hurst exponent
% Lx: length of topography in x direction.
% m: number of pixels in x
% n: number of pixels in y
% qr: roll-off wavelength
% qs: small wavelength cut-off
% qL: large wavelength cut-off


% =========================================================================
% Check number of inputs
if nargin < 6 || nargin > 8
    error(['The code requires at least 6, and maximum 8, inputs.',...
        'Check the code description']);
end
% =========================================================================
% make surface size an even number
if mod(n,2)
    n = n -1;
end
if mod(m,2)
    m = m -1;
end
% =========================================================================
PixelWidth = Lx/m;

Lx = m * PixelWidth; % image length in x direction
Ly = n * PixelWidth; % image length in y direction

% =========================================================================
% Wavevectors (note that q = (2*pi) / lambda; where lambda is wavelength.
% qx
qx = zeros(m,1);
for k = 0:m-1
    qx(k+1)=(2*pi/m)*(k);
end
qx = (unwrap((fftshift(qx))-2*pi))/PixelWidth;

% qy
qy = zeros(n,1);
for k = 0:n-1
    qy(k+1)=(2*pi/n)*(k);
end
qy = (unwrap((fftshift(qy))-2*pi))/PixelWidth;

% 2D matrix of q values
[qxx,qyy] = meshgrid(qx , qy);
[~,rho] = cart2pol(qxx,qyy);

% =========================================================================
% handle qr case
 if ~exist('qr','var')
     % default qr
      qr = 0; % no roll-off
 end
 if ~exist('qL','var')
     % default qL
      qL = 0; % no roll-off
 end
 qs = 2*pi*qs/Lx;
 qr = 2*pi*qr/Lx;
 qL = 2*pi*qL/Lx;
 
% 2D matrix of Cq values
Cq = zeros(n,m);
for i = 1:m
    for j = 1:n
        if rho(j,i) < qr
            Cq(j,i) = 0;
        elseif rho(j,i) < qL
            Cq(j,i) = qL^(-2*(H+1));
        elseif rho(j,i) < qs
            Cq(j,i) = rho(j,i).^(-2*(H+1));
        else
            Cq(j,i) = 0;
        end
    end
end

% =========================================================================
% applying rms
Cq(n/2+1,m/2+1) = 0; % remove mean
RMS_F2D = sqrt((sum(sum(Cq)))*(((2*pi)^2)/(Lx*Ly)));
alfa = sigma/RMS_F2D;
Cq = Cq.*alfa^2;

% =========================================================================
% Radial Averaging : DATA WILL BE USEFUL FOR CHECKING RESULTS
rhof = floor(rho);
J = 200; % resolution in q space (increase if you want)
qrmin = log10(sqrt((((2*pi)/Lx)^2+((2*pi)/Ly)^2)));
qrmax = log10(sqrt(qx(end).^2 + qy(end).^2)); % Nyquist
q = floor(10.^linspace(qrmin,qrmax,J));

C_AVE = zeros(1,length(q));
ind = cell(length(q)-1 , 1);
for j = 1:length(q)-1
    ind{j} = find( rhof > q(j) & rhof <=(q(j+1)));
    C_AVE(j) = nanmean(Cq(ind{j}));
end
ind = ~isnan(C_AVE);
C = C_AVE(ind);
q = q(ind);

% =========================================================================
% reversing opertation: PSD to fft
Bq = sqrt(Cq./(PixelWidth^2/((n*m)*((2*pi)^2))));

% =========================================================================
% apply conjugate symmetry to magnitude
Bq(1,1) = 0;
Bq(1,m/2+1) = 0;
Bq(n/2+1,m/2+1) = 0;
Bq(n/2+1,1) = 0;
Bq(2:end,2:m/2) = rot90(Bq(2:end,m/2+2:end),2);
Bq(1,2:m/2) = rot90(Bq(1,m/2+2:end),2);
Bq(n/2+2:end,1) = rot90(Bq(2:n/2,1),2);
Bq(n/2+2:end,m/2+1) = rot90(Bq(2:n/2,m/2+1),2);

% =========================================================================
% defining a random phase between -pi and phi (due to fftshift,
% otherwise between 0 and 2pi)
phi =  -pi + (pi+pi)*rand(n,m);

% =========================================================================
% apply conjugate symmetry to phase
phi(1,1) = 0;
phi(1,m/2+1) = 0;
phi(n/2+1,m/2+1) = 0;
phi(n/2+1,1) = 0;
phi(2:end,2:m/2) = -rot90(phi(2:end,m/2+2:end),2);
phi(1,2:m/2) = -rot90(phi(1,m/2+2:end),2);
phi(n/2+2:end,1) = -rot90(phi(2:n/2,1),2);
phi(n/2+2:end,m/2+1) = -rot90(phi(2:n/2,m/2+1),2);

% =========================================================================
% Generate topography
[a,b] = pol2cart(phi,Bq);
Hm = complex(a,b); % Complex Hm with 'Bq' as abosulte values and 'phi' as
% phase components.

z = ifft2(ifftshift((Hm))); % generate surface

PSD.qx = qx;
PSD.qy = qy;
PSD.Cq = Cq;
PSD.q = q;
PSD.C = C;


