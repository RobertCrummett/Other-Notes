% Plot IGRF model
clear; clf; close all; clc

% Load kmz data, kmz2struct
addpath(".\kmz2struct\DATA\","kmz2struct\")

% Open continent data
fn = "CONTDAT.kml";
[LAT,LON] = kmzLatLon(fn); clear fn

LAT = deg2rad(LAT); LON = deg2rad(LON);
[XL,YL,ZL] = sph2cart(LON,LAT,1.001);

% Load IGRF Coefficients
opts = spreadsheetImportOptions("NumVariables", 29);
opts.Sheet = "IGRF13coeffs";
opts.DataRange = "A5:AC199";
opts.VariableNames = ["cossin", "degree", "order", "IGRF", "IGRF1",...
    "IGRF2", "IGRF3", "IGRF4", "IGRF5", "IGRF6", "IGRF7", "IGRF8",...
    "DGRF", "DGRF1", "DGRF2", "DGRF3", "DGRF4", "DGRF5", "DGRF6",...
    "DGRF7", "DGRF8", "DGRF9", "DGRF10", "DGRF11", "DGRF12", "DGRF13",...
    "DGRF14", "IGRF9", "SV"];
opts.VariableTypes = ["categorical", "double", "double", "double",...
    "double", "double", "double", "double", "double", "double",...
    "double", "double", "double", "double", "double", "double",...
    "double", "double", "double", "double", "double", "double",...
    "double", "double", "double", "double", "double", "double", "double"];
opts = setvaropts(opts, "cossin", "EmptyFieldRule", "auto");
name = "IGRF13coeffs.xls";
IGRF13coeffs = readtable(name, opts, "UseExcel", false);
clear opts name

YEARKEY = dictionary(1900:5:2020,["IGRF", "IGRF1",...
    "IGRF2", "IGRF3", "IGRF4", "IGRF5", "IGRF6", "IGRF7", "IGRF8",...
    "DGRF", "DGRF1", "DGRF2", "DGRF3", "DGRF4", "DGRF5", "DGRF6",...
    "DGRF7", "DGRF8", "DGRF9", "DGRF10", "DGRF11", "DGRF12", "DGRF13",...
    "DGRF14", "IGRF9"]);

% Radius of Earth
R = 6.371e6; 

% Evaluation surface - at surface for my purposes
A = 6.371e6;

% Choose to include or exclude dipole term
plotDipoleField = true;

if plotDipoleField
    FD = 1;
elseif ~plotDipoleField
    FD = 2;
end

% Degrees and orders of Gauss coefficients
SC = IGRF13coeffs.cossin; % cos/sin, or g/h reference
DEG = IGRF13coeffs.degree; % n
ORD = IGRF13coeffs.order; % m
YEARSTEP = 0.25; % year increment

% Preallocation of frame data space
idf = 1;
FRAME = struct('cdata', cell(1, length(1900:YEARSTEP:2020)),...
    'colormap', cell(1, length(1900:YEARSTEP:2020)));

% Spline IGRF coefficients together
COEFF = zeros(195,length(1900:YEARSTEP:2020));
for IDC = 1:195 % rows
    COEFF(IDC,:) = interp1(1900:5:2020,IGRF13coeffs{IDC,4:28},...
        1900:YEARSTEP:2020);
end

IDC = 1;
for YEAR = 1900:YEARSTEP:2020

    LOWER = 5*floor(YEAR/5);
    
    % Get Gauss coefficients
    COE = COEFF(:,IDC);
    IDC = IDC + 1;

    % Positions
    step = 0.03; % radians, decrease to increase resolution
    phi = -pi:step:pi+step; % longitude
    the = 0:step:pi; % colatitude
    
    % Evaluate magnetic field
    % Preallocate space for radial magnetic field
    phi = phi(:)'; the = the(:);
    [PHI,THE] = meshgrid(phi,the);
    BR = zeros(size(PHI));
    
    for n = FD:13
    NTERM = (n+1)*(A/R)^(n+2);
    P = legendre(n,cos(the),'sch');
    for m = 0:n
        D = DEG == n;
        O = ORD == m;
        IDX = (D + O) == 2;
        C = COE(IDX);
        PM = P(m+1,:)';
    
        switch sum(IDX)
            case 2
                G = C(1); H = C(2);
                MTERM = (G*cos(m*PHI) + H*sin(m*PHI)).*PM;
            case 1
                G = C(1);
                MTERM = G*cos(m*PHI).*PM;
        end
        BR = BR + NTERM*MTERM;
    end
    end
    
    % Preallocate space for lateral magnetic field
    BP = zeros(size(PHI));
    
    for n = FD:13
    NTERM = (A/R)^(n+1);
    P = legendre(n,cos(the),'sch');
    for m = 0:n
        D = DEG == n;
        O = ORD == m;
        IDX = (D + O) == 2;
        C = COE(IDX);
        PM = P(m+1,:)';
    
        switch sum(IDX)
            case 2
                G = C(1); H = C(2);
                MTERM = (-G*m*sin(m*PHI) + H*m*cos(m*PHI)).*PM;
            case 1
                G = C(1);
                MTERM = -G*m*sin(m*PHI).*PM;
        end
        BP = BP + NTERM*MTERM;
    end
    end
    
    % Preallocate space for long magnetic field
    BT = zeros(size(PHI));
    
    for n = FD:13
    NTERM = (A/R)^(n+1);
    P = d_legendre(n,cos(the),'sch');
    for m = 0:n
        D = DEG == n;
        O = ORD == m;
        IDX = (D + O) == 2;
        C = COE(IDX);
        PM = P(m+1,:)';
    
        switch sum(IDX)
            case 2
                G = C(1); H = C(2);
                MTERM = (G*cos(m*PHI) + H*sin(m*PHI)).*PM;
            case 1
                G = C(1);
                MTERM = G*cos(m*PHI).*PM;
        end
        BT = BT + NTERM*MTERM;
    end
    end
    
    % TMI
    T = sqrt(BR.^2 + BT.^2 + BP.^2);
    
    BX = -BT;
    BY = BP;
    BZ = -BR;
    
    % Horizontal intensity
    H = sqrt(BX.^2 + BY.^2);
    % Inclination
    I = atan2(BZ, H);
    % Declination
    D = asin(BY./H);
    
    figure(1)
    set(gcf,'Color','white')
    CL = [2.28e4 6.91e4]; % colorbar limits, determined iteratively
    plotSphere(PHI,THE,T,CL)

    hold on

    % Draw points for colomap
    patch([0 0], [0 0], [min(T(:)) max(T(:))])
    colormap jet
    hCB = colorbar;
    hCB.Title.String='nT';
    hCB.Title.FontSize = 10;
    clim([CL(1) CL(2)])
    hCB.TickLength = 0;
    
    % Plot continent borders
    plot3(XL,YL,ZL,'k-')
    hold off

    view(-40,5) % POV
    axis equal; grid off; axis off
    title(strcat("Total Intensity, IGRF ",num2str(LOWER)),'FontSize',20,...
        'FontName','TimesNewRoman')
    set(gcf,"Position",[492.2000 49.8000 913.6000 828.8000]) % Screen size
    drawnow
    
    % Save frame data
    FRAME(idf) = getframe(gcf);
    idf = idf + 1;

end

% Write collection of frames to mp4
fn = "Temporary.mp4"; % File name / location
try
    warning('off'); delete(fn); warning('on')
catch
end
video_file = VideoWriter(fn, "MPEG-4");
video_file.FrameRate = 15; 
open(video_file)
writeVideo(video_file, FRAME)
close(video_file)