%% Tutorial for the third-party modified version QuaDRiGa-stat of the Fraunhofer QuaDRiGa
% Author: Lorenzo Miretti
%
%   This tutorial shows how to use QuaDRiGa-stat to generate i.i.d. channel
%   samples that can be used to perform conventional statistical
%   simulations of MIMO fading channels. In addition, it shows how obtain
%   relevant statistical parameters of the fading process such as the
%   spatial channel covariance matrix. The tutorial focuses on a
%   distributed MIMO setup with multi-antenna UEs and multi-antenna APs. 
%
%   Tested with MATLAB R2023b, 3GPP baseline mode, and single snapshot
%   layouts.

%% Initialization
clear variables
close all
addpath quadriga_src

%% Set up simulation parameters
s = qd_simulation_parameters;                           % Set up simulation parameters
s.show_progress_bars = 0;                               % Disable progress bars
s.center_frequency = 2.53e9;                            % Set center frequency
s.use_3GPP_baseline = 0;                                % Disable spherical waves and geometric polarization 
s.use_absolute_delays = 1;                              % Include delay of the LOS path. This is important for distributed MIMO setups.

%% Set up antennas
a_omni = qd_arrayant('omni');                                   % omnidirectional antenna

% AP antenna
n_vertical_AP = 4;
n_horizontal_AP = 4;
a_AP = qd_arrayant( '3gpp-3d', n_vertical_AP, ...               % UPA with vertically polarized antennas (pol_indicator = 1). 
    n_horizontal_AP, s.center_frequency, 1); 
a_AP.Fa = a_omni.Fa(:,:,ones(1,n_horizontal_AP*n_vertical_AP)); % Replace 3GPP antenna patterns with omnidirectional pattern.
a_AP.Fb = a_omni.Fb(:,:,ones(1,n_horizontal_AP*n_vertical_AP)); % This is to avoid that the array is shielded from the back.

% Pairs of cross polarized antennas can be obtained by setting pol_indicator = 3, 
% but the antenna patterns should be kept as they are, or replaced with
% qd_arrayant('xpol'), or similar.

% UE antenna
n_vertical_UE = 1;
n_horizontal_UE = 2;
a_UE = qd_arrayant( '3gpp-3d', n_vertical_UE, ...
    n_horizontal_UE, s.center_frequency, 1);
a_UE.Fa = a_omni.Fa(:,:,ones(1,n_horizontal_UE*n_vertical_UE));
a_UE.Fb = a_omni.Fb(:,:,ones(1,n_horizontal_UE*n_vertical_UE));


%% Set up layout

% By default, QuaDRiGa uses the downlink convention for naming and MIMO
% notation. This can be changed later.
l = qd_layout(s);                                       % Create new QuaDRiGa layout
l.set_scenario('BERLIN_UMa_NLOS');                      % Set common propagation scenario for all links (could be modified individually)
l.no_tx = 16;                                           % Set number of APs 
if not(sqrt(l.no_tx) == round(sqrt(l.no_tx)))           % Check if not a square number
    error('The number of APs must be a square number.')
end
l.no_rx = 64;                                            % Set number of UEs

% Set AP positions: uniform grid of access points in a squared service area
len = 500;                                              % Set lenght of squared service area [m]
x_grid = len/sqrt(l.no_tx)/2:len/sqrt(l.no_tx):len;     % Generate x coordinates of grid  
[x_grid,y_grid] = meshgrid(x_grid);                     % Generate coordinates of 2D grid 
l.tx_position(1:2,:) = [x_grid(:),y_grid(:)].';         % Set 2D AP coordinates  
l.tx_position(3,:) = 10;                                % Set height of APs [m]

% Set random UE positions: uniformly distributed in the above service area
l.rx_position(1:2,:) = len*rand(2,l.no_rx);             % Generate random 2D UE positions                                                 
l.rx_position(3,:) = 1.5;                               % Set height of UEs [m]

% Assign antennas and apply random rotations
for ii = 1:l.no_tx
    l.tx_array(ii) = copy(a_AP);                        % Assign copy of antenna to each AP
    l.tx_array(ii).rotate_pattern(360*(rand(1)),'z')    % Apply random azimuth rotation 
end
for ii = 1:l.no_rx
    l.rx_array(ii) = copy(a_UE);                        % Assign copy of antenna to each UE
    l.rx_array(ii).rotate_pattern(360*(rand(1)),'z')    % Apply random azimuth rotation 
    l.rx_array(ii).rotate_pattern(180*(rand(1)),'y')    % Apply random elevation rotation
end

% Plot layout
l.visualize([],[],0);                                   

% TODO For very large scale setups, modify pairing matrix (see tutorial 6) 
% to avoid computing and storing channels of very distant links
% The function set_pairing could be used to this end

%% Generate channels

% At the moment, the following function works only works for layouts with a
% single snapshot, as above. TODO (low priority): fix multi-snapshot case.
c = l.get_channels_subpaths();                            
% Remove per-antenna delays (3GPP mode does not generate them anyway) 
for tx = 1:l.no_tx
    for rx = 1:l.no_rx
        c(rx,tx).individual_delays = 0;                 
    end
end 
% Swap channels to focus on uplink convention  
swap_tx_rx(c);                                                                    

%% Example 1: narrowband or per-subcarrier simulations
ue = 1; ap = 1;                                                 % Focus on an arbitrary AP-UE pair

% Compute large scale fading parameters
beta = norm(c(ap,ue).coeff,"fro")^2;                            % Compute path loss
[K,K_AP,K_UE] = c(ap,ue).covariance();                          % Compute full, transmit, and receive spatial covariance matrices
% At the moment, only the zero mean case is supported 
% TODO: implement non-zero mean case  (fixed LoS phase)

% Draw an i.i.d. channel realization on an arbitrary subcarrier
randomize_phases(c(ap,ue));                                     % Draw i.i.d. phases for each path
H = sum(c(ap,ue).coeff,3);                                      % Sum over all paths 

%% Example 2: wideband simulations

% Set wideband system parameters
bandwidth = 10^8;                                               % 100 MHz bandwidth
n_subcarriers = 512;                                            % Number of subcarriers
subcarriers = linspace(-0.5,0.5,n_subcarriers);                 % Discrete frequency grid

% Compute large scale fading parameters
beta = norm(c(ap,ue).coeff,"fro")^2;                            % Compute path loss
[K,K_AP,K_UE] = c(ap,ue).covariance(bandwidth,subcarriers);     % Compute full, transmit, and receive space-frequency covariance matrices
% Note: only the first columns are computed and stored. 
% The remain parts can be recovered since the matrices are block Toeplitz
% and Hermitian symmetric.

% Draw an i.i.d. wideband channel realization 
randomize_phases(c(ap,ue));                                     % Draw i.i.d. phases for each path
H = c(ap,ue).fr(bandwidth,subcarriers);                         % Get frequency response 






