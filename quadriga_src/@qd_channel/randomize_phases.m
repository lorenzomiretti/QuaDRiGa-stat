function randomize_phases( h_channel )
%RANDOMIZE_PHASES applies an i.i.d. phase shift to all paths.
%
% Calling object:
%   Single object
%
% Description:
%   This method updates and returns the channel coefficients by applying an
%   i.i.d. phase shift to all paths. If the channel has been generated
%   using the get_channel_subpaths() method, then this method can be
%   used to draw a new small-scale fading realization caused by microscopic
%   mobility (on the order of the wavelenght) following a statistical 
%   approach. 
%   
%   Note: all snapshots are updated using the same phase shift, to keep the
%   deterministic phase relations across multiple snapshots consistent. 



if numel( h_channel ) > 1
    error('QuaDRiGa:qd_channel:randomize_phases','??? "randomize_phases" is only defined for scalar objects.')
else
    h_channel = h_channel(1,1); % workaround for octave
end

% Get the channel tensor
H = h_channel.coeff;
phi = 2*pi*rand(h_channel.no_path);
exp_jphi = exp(1j*phi);
% Update channel tensor
for i_snap = 1:h_channel.no_snap
    for i_path = 1:h_channel.no_path
        H(:,:,i_path,i_snap) = H(:,:,i_path,i_snap)*exp_jphi(i_path); 
    end
end
h_channel.coeff = H;
end
