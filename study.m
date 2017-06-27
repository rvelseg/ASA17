%% This code is distributed under the terms of the Berkeley Software
%% Distribution (BSD) license
%% http://www.opensource.org/licenses/bsd-license.php

%% Copyright (c) 2017, Roberto Velasco Segura and Pablo L. Rendón
%% All rights reserved.

%% Redistribution and use in source and binary forms, with or without
%% modification, are permitted provided that the following conditions are
%% met:

%% 1. Redistributions of source code must retain the above copyright
%% notice, this list of conditions and the following disclaimer.

%% 2. Redistributions in binary form must reproduce the above copyright
%% notice, this list of conditions and the following disclaimer in the
%% documentation and/or other materials provided with the distribution.

%% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
%% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
%% LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
%% A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
%% HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
%% SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
%% LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
%% DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
%% THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
%% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
%% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function study(varargin)

  %% Exceution of this script took about 10 hours.

  %% Experiment number
  n = [1:43];

  if length(varargin) >= 1
    i_ini = varargin{1};
  else
    i_ini = length(n);
  end

  %% 10 --> 3200
  k = round(10.*(1.15.^(n-1)));  %% 1.15 ~ 2**(1/5)

  %% coeficient to amplify feedback
  A = linspace(10,50,length(n));

  t_max(1:17)         = 10;
  t_max(18:length(n)) = 10;

  %% There is acctually another level of the study, where you vary the
  %% following parameters, by hand. A fully automated study would have
  %% arrays here, and iterate over them. In that case, something
  %% should be added 'output_file_prefix' to distinguish results.
  NN = 1;
  f_op = 2;  %% only values 2 or 4 are allowed
  fixed_width = true;

  for i=i_ini:-1:1

    close all;

    i_str=sprintf('%.3d',i);

    %% A very rough estimation of the oscilation period. It is
    %% important for this value to be larger than the real period.
    tau = 50/n(i);

    c = 1;
    x_max = 1;

    Nx = 100;
    CFL = 0.9;
    dx = x_max/(Nx - 1);
    dt = CFL*dx/c;

    mic_pos = [4*dx];                  % left-most point, not in the bc
    mic_pos = [mic_pos, 2/3];          % antinodes of the second mode
    mic_pos = [mic_pos, 2/5, 4/5];     % antinodes of the third mode
    mic_pos = [mic_pos, 2/7, 4/7, 6/7];
    mic_pos = [mic_pos, 2/9, 4/9, 6/9, 8/9];

    props = struct( ...
		  'ampl', 1, ...         % amplitude of the emitted wave
		  'x_max', x_max, ...    % domain length in meters
		  'Nx', Nx, ...          % number of grid points
		  't_max', t_max(i), ... % final time in seconds
		  'c', c, ...            % speed of sound in m/s
		  'CFL', CFL, ...        % Courant-Friedrichs-Lewy constant
		  'delta', [0.25, 0.35], ...
		  'att_scheme', [1, 2], ...
		  'split', true, ...
		  'mic_pos', mic_pos, ...
		  'mic_plot_png', false, ...
		  'display_every', round((tau/dt)/20)+1, ...
		  'display_ylim', 15, ...
		  'snaps_t', linspace(t_max(i)-tau, t_max(i), 20), ...
		  'output_dir', 'results', ...
		  'output_file_prefix', ['lips_', i_str, '_'], ...
		  'video_framerate', 5, ...
		  'video_file', true ...
		  );

    %% The parameter named 'something' is there just for testing the
    %% structure.
    props.study_params = struct( ...
				 'n', n(i), ...
				 'k', k(i), ...
				 'tau', tau, ...
				 'A', A(i), ...
				 'NN', NN, ...
				 'f_op', f_op, ...
				 'fixed_width', fixed_width ...
			       );

    disp('------------------------------------');
    disp(['Performing simulation ', num2str(i)])
    disp(props)
    disp(props.study_params)

    simulation(props);
    postprocessing1(props);
  end

  postprocessing2(props, n, k);
end
