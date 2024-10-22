function [K, K_tx, K_rx] = covariance( h_channel, bandwidth, carriers)
%COVARIANCE returns the channel covariance matrices
%
% Calling object:
%   Single object
%
% Description:
%   This method updates and returns the channel covariance matrices
%   assuming that new channel realizations are drawn using the
%   randomized_phases method, that channels are generated using the
%   get_channels_subpaths method, and identical delays across the MIMO
%   array, i.e., individual_delays = 0. If bandwidth and carriers are
%   omitted, the method returns the spatial covariance matrix. Otherwise,
%   it returns the space-frequency covariance matrices.
%   
% Input:
%   bandwidth
%       Hz
%   carriers
%       Row vector with values in [-1/2, 1/2)
% Output:
%   K, K_tx, K_rx
%       The full, transmit, and receive channel covariance matrices. In the
%       spatial only case, the matrices have dimension [no_ant, no_ant],
%       where no_antennas equals no_txant*no_rxant, no_txant, no_rxant,
%       respectively. In the space-frequency case, the matrices have
%       dimension [no_ant*no_subcarriers, no_ant]. The rest of the covariance
%       matrices can be recovered since they are block Toeplitz Hermitian symmetric. 


if numel( h_channel ) > 1
    error('QuaDRiGa:qd_channel:covariance','??? "covariance" is only defined for scalar objects.')
else
    h_channel = h_channel(1,1); % workaround for octave
end

if ~exist( 'bandwidth','var' ) || isempty( bandwidth ) || ...
        ~exist( 'carriers','var' ) || isempty(carriers)
    wideband = false;
else
    wideband = true;
end

% TODO check input as in the fr method

% Get the channel tensor
H = h_channel.coeff;
n_rx = h_channel.no_rxant;
n_tx = h_channel.no_txant;
n_path = h_channel.no_path;
n_snap = h_channel.no_snap;


if ~wideband
    % Compute spatial covariance matrices 
    K_tx = zeros(n_tx,n_tx,n_snap);                        % tx covariance matrix
    K_rx = zeros(n_rx,n_rx,n_snap);                        % rx covariance matrix
    K = zeros(n_rx*n_tx,n_rx*n_tx,n_snap);                 % full covariance matrix
    for i_snap = 1:n_snap
        for i_path = 1:n_path
            H_i = H(:,:,i_path,i_snap);
            K_rx(:,:,i_snap) = K_rx(:,:,i_snap) + H_i*H_i';
            K_tx(:,:,i_snap) = K_tx(:,:,i_snap) + H_i'*H_i;
            K(:,:,i_snap) = K(:,:,i_snap) + H_i(:)*H_i(:)';
        end
    end
else
    n_carrier = length(carriers);
    % Get the delays
    tau = h_channel.delay;
    % Normalize delays
    d = tau.*bandwidth;

    % Compute space-frequency covariance matrices 
    % For efficiency resons, we keep only the first n_ant columns, since 
    % the matrices are block Toeplitz Hermitian symmetric
    K_tx = zeros(n_tx*n_carrier,n_tx,n_snap);                        % tx covariance matrix
    K_rx = zeros(n_rx*n_carrier,n_rx,n_snap);                        % rx covariance matrix
    K = zeros(n_rx*n_tx*n_carrier,n_rx*n_tx,n_snap);                 % full covariance matrix
    for i_snap = 1:n_snap
        for i_path = 1:n_path
            % Compute spatial parts
            H_i = H(:,:,i_path,i_snap);
            G_i_rx = H_i*H_i';
            G_i_tx = H_i'*H_i;
            G_i = H_i(:)*H_i(:)';
            % Compute frequency correlation matrix
            T1 = exp(-2j*pi*(carriers-carriers(1))*d(i_path,i_snap)).';
            % Update covariance matrices
            K_tx(:,:,i_snap) = K_tx(:,:,i_snap)+kron(T1,G_i_tx);
            K_rx(:,:,i_snap) = K_rx(:,:,i_snap)+kron(T1,G_i_rx);
            K(:,:,i_snap) = K(:,:,i_snap)+kron(T1,G_i);  
        end
    end  
end
