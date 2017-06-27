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

function postprocessing1(props)

  file = [props.output_dir, ...
	  filesep, ...
	  props.output_file_prefix, ...
	  'mic.csv'];
  mic = csvread(file);

  CFL = props.CFL;
  Nx = props.Nx;
  x_max = props.x_max;
  c = props.c;
  dx = x_max/(Nx-1.0);
  dt = CFL*dx/c;

  mic_amount = length(props.mic_pos);

  t_1 = {};
  y_1 = {};
  t_2 = {};
  y_2 = {};
  f_1 = {};
  Y_1 = {};
  f_2 = {};
  Y_2 = {};
  f_pks = {};
  Y_pks = {};

  A_a = [];
  %% pks_min_a_diff = [];
  pks_a = {};
  pks_a_i = {};
  Y_a = [];
  f_a = [];

  for i=1:mic_amount

    y_o = mic(:,i);

    t_o = [0.0:dt:dt*(length(y_o)-1)];
    %% length of t must be equal to the length of y
    if length(t_o) ~= length(y_o)
      disp('Unconsistent data')
      return
    end

    %% select a susbset of the time series, first half
    t_max = max(t_o);
    t_s = find(t_o > 0 & t_o < t_max/2);
    t = t_o(t_s);
    y = y_o(t_s);

    %% ensure an even number of measures
    if mod(length(y),2) == 1
      y = y(1:end-1);
      t = t(1:end-1);
    end
    N = length(y);

    t_1{i} = t_o;
    y_1{i} = y_o;

    t_p_i = find(t_o > t_max - props.study_params.tau);
    t_2{i} = t_o(t_p_i);
    y_2{i} = y_o(t_p_i);

    f = (1.0/dt)*(0:(N/2))/N;

    Y = fft(y)*(2/N);
    Y = Y(1:N/2+1);

    f_1{i} = f;
    Y_1{i} = abs(Y);

    %% select a susbset of the time series, second half
    t_s = find(t_o > (t_max/2) & t_o < t_max);
    t = t_o(t_s);
    y = y_o(t_s);
    A_a(i) = max(y) - min(y);

    %% ensure an even number of measures
    if mod(length(y),2) == 1
      y = y(1:end-1);
      t = t(1:end-1);
    end
    N = length(y);

    f = (1.0/dt)*(0:(N/2))/N;

    Y = fft(y)*(2/N);
    Y = Y(1:N/2+1);

    f_2{i} = f;
    Y_2{i} = abs(Y);

    %% Peaks
    is_octave = ( exist('OCTAVE_VERSION', 'builtin') > 0 );
    if is_octave
      try
        pkg load signal;
      catch E
    	disp(['ERROR: You need to install and/or', ...
	      'load the "signal" package for octave.'])
	rethrow(E)
      end
    end

    df = f(2)-f(1);
    [pks, pks_i] = findpeaks(abs(Y), ...
			     'MinPeakDistance', (0.5/df)*0.6);

    f_pks{i} = f(pks_i);
    Y_pks{i} = pks;

    %% [pks_max, pks_max_i] = max(pks);
    %% pks_min_diff = pks_max;
    %% %% pks_min_diff_i = NaN;
    %% for ip=1:length(pks)
    %%   if pks_min_diff > abs(pks_max - pks(ip)) && ...
    %% 	 pks_max_i ~= ip
    %% 	pks_min_diff = abs(pks_max - pks(ip));
    %% 	%% pks_min_diff_i = pks_i(ip);
    %%   end
    %% end
    %% pks_min_a_diff(i) = pks_min_diff;

    pks_a{i} = pks;
    pks_a_i{i} = pks_i;
    Y_a(i,:) = Y;
    f_a(i,:) = f;

    figure_file_name =  [props.output_dir, ...
			 filesep, ...
			 props.output_file_prefix, ...
			 'signal_n_fft_mic_', ...
			 num2str(i), ...
			 '.png'];

    make_plot(t_1{i},   y_1{i}, ...
	      t_2{i},   y_2{i}, ...
	      f_1{i},   Y_1{i}, ...
	      f_2{i},   Y_2{i}, ...
	      f_pks{i}, Y_pks{i}, ...
	      t_max, i, props, figure_file_name);
  end

  %% Choose the best microphone, but not the first one, since this one
  %% is to close to the source.
  [pks_best_a_v, pks_best_a_i] = max(A_a(2:end));
  pks_best_a_i = pks_best_a_i + 1;  % add 1 to get the index over A_a

  pks_best   = pks_a{pks_best_a_i};
  pks_best_i = pks_a_i{pks_best_a_i};
  Y_best = Y_a(pks_best_a_i,:);
  f_best = f_a(pks_best_a_i,:);

  csvwrite([props.output_dir, ...
  	    filesep, ...
  	    props.output_file_prefix, ...
  	    'peaks_best.csv'], ...
  	   [pks_best_i, pks_best, pks_best_a_i*ones(length(pks_best),1)]);
  csvwrite([props.output_dir, ...
  	    filesep, ...
  	    props.output_file_prefix, ...
  	    'f_Y.csv'], ...
  	   [f_best', Y_best']);

  figure_file_name =  [props.output_dir, ...
		       filesep, ...
		       props.output_file_prefix, ...
		       'signal_n_fft_mic_best', ...
		       '.png'];

  make_plot(t_1{pks_best_a_i},   y_1{pks_best_a_i}, ...
	    t_2{pks_best_a_i},   y_2{pks_best_a_i}, ...
	    f_1{pks_best_a_i},   Y_1{pks_best_a_i}, ...
	    f_2{pks_best_a_i},   Y_2{pks_best_a_i}, ...
	    f_pks{pks_best_a_i}, Y_pks{pks_best_a_i}, ...
	    t_max, pks_best_a_i, props, figure_file_name);
end


function make_plot(t_1, y_1, ...
		   t_2, y_2, ...
		   f_1, Y_1, ...
		   f_2, Y_2, ...
		   f_pks, Y_pks, ...
		   t_max, i, props, figure_file_name)
  figure(50)
  clf(50)

  subplot(3,2,1);
  plot(t_1,y_1)
  xlabel('time (s)')
  ylabel('y')
  xlim([0, t_max])
  title(['mic ', num2str(i), ...
	 ', on position ', num2str(props.mic_pos(i)), ...
	 ', n = ', num2str(props.study_params.n)])

  subplot(3,2,2);
  plot(t_2,y_2)
  xlabel('time (s)')
  ylabel('y')
  xlim([t_max - props.study_params.tau, t_max])

  subplot(3,1,2);
  stem(f_1,Y_1)
  title('Spectrum, first half of the timeseries.')
  xlabel('Frequency (Hz)')
  ylabel('|Y|')
  xlim([0,6])
  set(gca,'Xtick',0:0.25:6)
  grid on;
  set(gca,'box','on')
  %% set(gca,'Yscale','log')

  subplot(3,1,3);
  hold on;
  stem(f_2,Y_2)
  plot(f_pks, Y_pks, 'v', 'MarkerSize',12)
  hold off;
  title('Spectrum, second half of the timeseries.')
  xlabel('Frequency (Hz)')
  ylabel('|Y|')
  xlim([0,6])
  grid on;
  set(gca,'box','on')
  set(gca,'Xtick',0:0.25:6)
  %% set(gca,'Yscale','log')

  print(figure_file_name, '-dpng');
end
