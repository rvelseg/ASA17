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

function postprocessing2(props, n, k)

  input1 = ['./results'];

  f = {};
  Y = {};
  f_i_m = NaN;

  figure(5001)

  for i=1:length(n)

    i_str=sprintf('%.3d',i);
    output_file_prefix = ['lips_', i_str, '_'];

    file = [input1, ...
	    filesep, ...
	    output_file_prefix, ...
	    'f_Y.csv'];

    data = csvread(file);
    f_line = data(:,1);
    Y_line = data(:,2);

    f{i} = f_line;
    Y{i} = Y_line;

    if isnan(f_i_m)
      f_i_m = length(find(f_line < 4));
    else
      f_i_m = max(f_i_m, length(find(f_line < 4)));
    end
  end

  for i=1:length(n)
    Y_log = log(eps + abs(Y{i}(1:f_i_m)));
    xx(:,i) = k(i)*ones(1,f_i_m);
    yy(:,i) = f{i}(1:f_i_m);
    zz(:,i) = Y_log(1:f_i_m);
    zz2(:,i) = abs(Y{i}(1:f_i_m));
    zz3(:,i) = abs(Y{i}(1:f_i_m))/max(abs(Y{i}(1:f_i_m)));
  end

  clf(5001)
  surf(xx,yy,zz, ...
       'EdgeColor','none', ...
       'LineStyle','none', ...
       'FaceLighting','phong');
  xlabel('Stiffness s')
  ylabel('freq')
  view(2)
  set(gca, 'XScale', 'log')
  %% set(gca, 'XScale', 'linear')
  ylim([0,4])
  xlim([k(1), k(end)]);
  shading interp
  colormap jet
  figure_file_name = [props.output_dir, ...
  		      filesep, ...
  		      'k_f_log.png'];
  print(figure_file_name, '-dpng');

  clf(5001)
  surf(xx,yy,zz2, ...
       'EdgeColor','none', ...
       'LineStyle','none', ...
       'FaceLighting','phong');
  xlabel('Stiffness s')
  ylabel('freq')
  view(2)
  set(gca, 'XScale', 'log')
  %% set(gca, 'XScale', 'linear')
  ylim([0,4])
  xlim([k(1), k(end)]);
  shading interp
  colormap gray
  colormap(flipud(colormap))
  figure_file_name = [props.output_dir, ...
  		      filesep, ...
  		      'k_f.png'];
  print(figure_file_name, '-dpng');

  clf(5001)
  surf(xx,yy,zz3, ...
       'EdgeColor','none', ...
       'LineStyle','none', ...
       'FaceLighting','none');
  xlabel('Stiffness s')
  ylabel('freq')
  view(2)
  set(gca, 'XScale', 'log')
  %% set(gca, 'XScale', 'linear')
  ylim([0,4])
  xlim([k(1), k(end)]);
  shading interp
  colormap gray
  colormap(flipud(colormap))
  figure_file_name = [props.output_dir, ...
		      filesep, ...
		      'k_f_norm.png'];
  print(figure_file_name, '-dpng');

end
