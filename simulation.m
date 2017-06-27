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

function simulation(props)

  app = propagator(props);
  app.run(@app_ic, @app_bc);
end

%% Initial condition
function [y_now y_prev] = app_ic(obj)

  y_now = obj.y_now;
  y_prev = obj.y_prev;

  for i = [1:length(obj.x)]
    y_now(i) = 0;
  end

  for i = [1:length(obj.x)]
    y_prev(i) = 0;
  end
end

%% Boundary condition
function [y_prev, y_now, y_next] = app_bc(obj)

  y_next = obj.y_next;
  y_now  = obj.y_now;
  y_prev = obj.y_prev;

  %% harmonic oscillator constants
  k = obj.study_params.k;           % Stiffness
  m = 8;             % Mass
  R = 0;

  %% Frequency of the undamped harmonic oscillator. This value is not
  %% used as you can see.
  freq_0 = sqrt(k/m);
  %% msg(obj,['freq_0 = ', num2str(freq_0)]);

  %% reflective boundary, Dirichlet
  y_next(end) = 0;

  valve_eq = 0;
  valve_max = 0.1;
  S = 1;

  %% The variable 'h' in the poster here is called valve_
  persistent valve_next;
  persistent valve_now;
  persistent valve_prev;
  persistent valve_start;
  persistent valve_period_m;
  persistent valve_period_m_f;
  if obj.first
    %%    disp('Initializing persistent variables');
    valve_next  = valve_eq-valve_max/2;
    valve_now   = valve_eq-valve_max/2;
    valve_prev  = valve_eq-valve_max/2;
    valve_start = 0;
    valve_period_m = 0;
    valve_period_m_f = 0;
  end

  valve_refinement = 50;
  for iv=1:valve_refinement
    valve_dt = obj.dt/valve_refinement;
    valve_next = ( + valve_now * ( + (R/valve_dt)  ...
				   + (2*m/(valve_dt^2)) ) ...
		   - valve_prev * ( m/(valve_dt^2) ) ...
		   + valve_eq * k ...
		   + y_next(1)^obj.study_params.NN ...  % n in the poster
		 ) / ...
		 ( + k ...
		   + (R/valve_dt) ...
		   + (m/(valve_dt^2)) );

    %% Do not remove this block, you will get (almost exclusivelly)
    %% solutions with unlimited exponential growth. Moreover,
    %% eventually the code crashes.
    if valve_next - valve_eq > valve_max
      valve_next = valve_eq + valve_max;
    elseif valve_next - valve_eq < -valve_max
      valve_next = valve_eq - valve_max;
    end

    %% if the line between valve_next and valve_prev cross the x axis
    if valve_next*valve_prev < 0
      valve_start = 1;
      if valve_period_m_f == 0
	valve_period_m_f = 1;
      elseif valve_period_m_f == 1
	valve_period_m_f = 2;
	valve_period_m = obj.n;
      elseif valve_period_m_f == 2
	valve_period_m_f = 3;
	valve_period_m = obj.n - valve_period_m;
      end
    end

    %% disp(['valve_next = ', num2str(valve_next)]);
    %% disp(['valve_now = ',  num2str(valve_now)]);
    %% disp(['valve_prev = ', num2str(valve_prev)]);

    valve_prev = valve_now;
    valve_now  = valve_next;
  end

  %% reflective boundary, Neumann
  y_now(1)  = y_now(2);
  y_next(1) = y_now(2);

  persistent mu
  persistent sigma

  if valve_start && valve_period_m_f == 3 %&& valve_next > 0

    %% Added p. See that after saturation option 1 and 2 are square
    %% pulses, but not option 3, or 4.
    if obj.study_params.f_op == 1
      disp('WARNING: option 1 for f(t) is not implemented.')
      %% %% Option 1, square root
      %% in = 0.3*obj.study_params.A ...
      %%      *(valve_next/abs(valve_next))...
      %%      *(abs(valve_next)^(0.5));
    elseif obj.study_params.f_op == 2
      %% Option 2, linear, square pulses
      in = obj.study_params.A*valve_next;
    elseif obj.study_params.f_op == 3
      disp('WARNING: option 3 for f(t) is not implemented.')
      %% %% Option 3, Gaussian pulses
      %% if valve_start == 1
      %%   mu = 0.4*valve_period_m;
      %%   sigma = valve_period_m*0.2;
      %% end
      %% g = 2*exp(-0.5*((valve_start - mu)/sigma)^2);
      %% in = obj.study_params.A*valve_next*g;
    elseif obj.study_params.f_op == 4
      %% Option 4, parabolic pulses
      if valve_start == 1
	mu = 0.4*valve_period_m;
      end
      pp = 2*(1-((mu-valve_start)/mu)^2);
      in = obj.study_params.A*valve_next*pp;
    end

    y_now(1) = y_now(1) + in;
    y_now(2) = y_now(2) + in;
    y_prev(1) = y_prev(1) + in;
    y_prev(2) = y_prev(2) + in;
    if obj.study_params.fixed_width
      valve_start = valve_start + 1;
      if valve_start > valve_period_m * 0.8
	valve_start = 0;
      end
    end
  end

  plot_valve(obj, ...
	     y_next, ...
	     valve_next, valve_eq, valve_max)
end

function msg(obj, message)
  if mod(obj.n,obj.display_every) == 0
    disp(message)
  end
end

function plot_valve(obj, ...
		    y_next, ...
		    valve_next, valve_eq, valve_max)

  persistent valve_v;
  persistent valve_hist;
  persistent y_diff_h;
  persistent hist_i;
  if obj.first
    if obj.video_file
      valve_v = VideoWriter([obj.output_dir, ...
			     obj.output_file_prefix, ...
			     'valve']);
      open(valve_v);
    end
    hist_i = 1;
    valve_hist = [];
    y_diff_h = [];
  end
  valve_hist(hist_i) = valve_next;
  y_diff_h(hist_i) = y_next(1);

  if mod(obj.n,obj.display_every) == 0
    figure(1101)
    clf(1101)
    hold on;
    %% plot([0,0], [0,valve_next], 'o');
    plot([1:hist_i], y_diff_h/100, 'b-', 'LineWidth', 2);
    plot([1:hist_i], valve_hist, 'r-', 'LineWidth', 2);
    hold off;
    set(gca, 'ylim', [valve_eq-valve_max*3, valve_eq+valve_max*3])
    legend('y/100', 'valve');
    grid on;

    if obj.video_file
      frame = getframe(1101);
      writeVideo(valve_v,frame);
    end
  end

  %% posn = get(gca, 'OuterPosition');
  %% posn(2) = posn(2) + 200;
  %% set(gca, 'OuterPosition', posn);

  hist_i = hist_i + 1;

  if obj.last
    if obj.video_file
      close(valve_v);
    end
  end
end
