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

classdef propagator

  properties

    %% ------- configurable parameters

    ampl = 1;     % amplitude of the emitted wave
    x_max = 1;    % domain length [m]
    Nx = 100;     % number of grid points
    t_max = 1;    % final time [s]
    c = 1;        % speed of sound in [m/s]
    CFL = 0.9;    % Courant-Friedrichs-Lewy constant
    delta = 0;    % Sound diffusivity coefficient

    %% ------- calculated and fixed paramenters

    dx = 0;          % spatial step [m]
    dt = 0;          % temporal step [s]
    t_min = 0.0;     % initial time [s]

    x = [];          % spatial domain
    t_axis = [];     % temporal domain

    %% ---------------- Cached variables

    cdtdx2 = 0;
    dt2 = 0;

    %% ---------------- Empty data structures for the soution

    n = 0;               % current step
    t = 0;               % time variable
    y_next = [];         % solution at t
    y_now  = [];         % solution at t-dt
    y_prev = [];         % solution at t-2*dt
    D = [];              % derivative dy/dt at t
    y_left = [];         % left going wave
    y_right = [];        % right going wave

    %% ----------------- Microphones

    mic_pos = [];
    mic = [];
    mic_num = 0;
    mic_pos_i = [];
    mic_csvfile = 'mic.csv';
    mic_plot_png = true;

    %% ----------------- Snapshots

    snaps = [];
    snaps_t = [];
    snaps_n = [];

    %% ----------------- Other internal variables

    first = 1;
    last  = 0;
    get_D = false;
    display_D = false;
    display_ylim = 2.2;
    display_D_ylim = NaN;
    display_every = 1;
    split = false;
    is_octave = ( exist('OCTAVE_VERSION', 'builtin') > 0 );
    video_file = false;
    video_framerate = NaN;
    matlabVerFlag1 = false; % use animatedline
    isok = true;
    att_scheme = 2;
    study_params = NaN;
    output_file_prefix = '';
    output_dir = '';

    %% collection of animatedline objects for main variables
    %% 1    y_now
    %% 2    y_left
    %% 3    y_right
    %% 4    D
    al = {};
    %% attributes for each line
    al_p = {{'LineStyle','-', ...   % 1
             'LineWidth',2, ...
             'Color','k'}, ...
            {'LineStyle','--', ...  % 2
             'LineWidth',2, ...
             'Color','g'}, ...
            {'LineStyle',':', ...  % 3
             'LineWidth',2, ...
             'Color','g'}};
  end

  methods
    function obj = propagator(varargin)

      if isstruct(varargin{1})
        props = varargin{1};
      elseif mod(length(varargin),2) == 0
        props = struct(varargin{:});
      else
        disp('Unrecognized options')
        varargin
        obj.isok = false;
        return
      end

      %% Get independent parameters, or use defalt values.
      if isfield(props,'ampl')
        obj.ampl = props.ampl;
      end
      if isfield(props,'x_max')
        obj.x_max = props.x_max;
      end
      if isfield(props,'Nx')
        obj.Nx = props.Nx;
      end
      if isfield(props,'t_max')
        obj.t_max = props.t_max;
      end
      if isfield(props,'c')
        obj.c = props.c;
      end
      if isfield(props,'CFL')
        obj.CFL = props.CFL;
      end

      %% calculate dependent parameters
      obj.dx = obj.x_max/(obj.Nx - 1.0);
      obj.dt = obj.CFL*obj.dx/obj.c;
      obj.x = linspace(0,obj.x_max,obj.Nx);  % spatial domain
      obj.cdtdx2 = (obj.c*obj.dt/obj.dx)^2;
      obj.dt2 = obj.dt^2;
      obj.t_axis = linspace(obj.t_min, ...
                            obj.t_max, ...
                            (obj.t_max-obj.t_min)/obj.dt);
      obj.y_next = zeros(obj.Nx,1);           % solution at t
      obj.y_now = obj.y_next;                 % solution at t-dt
      obj.y_prev = obj.y_now;                 % solution at t-2*dt

      %% get some other optional parameters
      if isfield(props, 'output_dir')
        l = length(props.output_dir);
        while ~isstrprop(props.output_dir(l),'alphanum')
          props.output_dir = props.output_dir(1:l-1);
          l = l - 1;
        end
        props.output_dir = [props.output_dir, filesep];
        obj.output_dir = props.output_dir;
        if ~exist(obj.output_dir,'dir')
          mkdir(obj.output_dir)
        end
      end
      if isfield(props, 'video_framerate')
        obj.video_framerate = props.video_framerate;
      end
      if isfield(props, 'output_file_prefix')
        obj.output_file_prefix = props.output_file_prefix;
      end
      if isfield(props, 'study_params')
        obj.study_params = props.study_params;
      end
      if isfield(props,'att_scheme')
        obj.att_scheme = props.att_scheme;
        for ai=1:length(obj.att_scheme)
          if ~ismember(obj.att_scheme(ai), [1,2,3])
            disp('ERROR: Unrecognized attenuation scheme: ')
            disp(obj.att_scheme(ai))
          end
        end
      end
      if isfield(props,'display_ylim')
        obj.display_ylim = props.display_ylim;
      else
        obj.display_ylim = 2.2*obj.ampl;
      end
      if isfield(props,'display_D')
        obj.display_D = props.display_D;
      end
      if isfield(props,'display_D_ylim')
        obj.display_D_ylim = props.display_D_ylim;
      end
      if isfield(props,'delta')
        obj.delta = props.delta;
      end
      if isfield(props,'display_every')
        obj.display_every = props.display_every;
      end
      if isfield(props,'mic_pos')
        obj.mic_pos = props.mic_pos;
        obj.mic_num = length(obj.mic_pos);
        obj.mic_pos_i = min(obj.Nx, ...
                            max(1, ...
                                round(obj.mic_pos*obj.Nx ...
                                      / obj.x_max)));
      end
      if isfield(props,'mic_csvfile')
        obj.mic_csvfile = props.mic_csvfile;
      end
      if isfield(props,'mic_plot_png')
        obj.mic_plot_png = props.mic_plot_png;
      end
      if isfield(props,'snaps_t')
        obj.snaps_t = props.snaps_t;
        obj.snaps_n = round(obj.snaps_t/obj.dt);
      end
      if isfield(props,'split')
        obj.split = props.split;
      end
      if obj.display_D || obj.split
        obj.get_D = true;
      end
      if isfield(props,'video_file')
        obj.video_file = props.video_file;
      end
      if obj.is_octave
        if obj.video_file
          disp(['Video file generation in octave ', ...
                'is not yet implemented, for this class.',...
                'The option "video_file" ', ...
                'is ignored.']);
          obj.video_file = false;
        end
      else
        obj.matlabVerFlag1 = ~verLessThan('matlab', '8.4');
      end
    end

    function run(obj, ic, bc)

      if ~obj.isok
        return
      end

      [obj.y_now obj.y_prev] = ic(obj);

      %% I don't really like these variables "figures_", an array of
      %% actual figures would be beter.
      figures = {};
      figures_ylim = {};
      figures_title = {};

      figures{1} = 1; % y_now plot
      figures_ylim{1} = [-obj.display_ylim, obj.display_ylim];
      figures_title{1} = 'Displacement y';

      if obj.split
        figures{1} = [figures{1}, 2, 3];  % y_left y_right plot
      end

      if obj.display_D
        figures{2} = 4;  % D plot
        if isnan(obj.display_D_ylim)
          obj.display_D_ylim = 2.2*max(obj.D_get());
        end
        figures_ylim{2} = [-obj.display_D_ylim, obj.display_D_ylim];
        figures_title{2} = 'Velocity dy/dt';
        disp('WARNING: figures could be overlapped.')
      end

      for i=1:length(figures)
        obj = obj.ani_init(i, ...
                           figures_ylim{i}, ...
                           figures{i}, ...
                           figures_title{i});
      end

      mic_i = 1;
      snaps_i = 1;

      if obj.video_file
        if ischar(obj.video_file)
          v = VideoWriter([obj.output_dir, ...
                           obj.output_file_prefix, ...
                           obj.video_file]);
        else
          v = VideoWriter([obj.output_dir, ...
                           obj.output_file_prefix, ...
                           'animation']);
        end
        if ~isnan(obj.video_framerate)
          v.FrameRate = obj.video_framerate;
        end
        open(v);
      end

      for t = obj.t_axis

        if t == obj.t_axis(end)
          obj.last = 1;
        end

        obj.t = t;
        if mod(obj.n,obj.display_every) == 0
          figures_y{1} = [];
          if obj.get_D
            obj.D = obj.D_get();
          end
          if obj.display_D
            obj.ani_step(2, ...
		         obj.D, ...
		         figures_ylim{2}, ...
		         4, ...
                         figures_title{2});
          end
          if obj.split
            [obj.y_left, obj.y_right] = obj.get_split();
            figures_y{1} = [figures_y{1}, obj.y_left, obj.y_right];
          end
          obj.ani_step(1, ...
		       [obj.y_now, figures_y{1}], ...
		       figures_ylim{1}, ...
		       figures{1}, ...
                       figures_title{1});

          if obj.video_file
            frame = getframe(1);
            writeVideo(v,frame);
          end
          obj.disp_time(t);
          drawnow;
          if obj.first
            pause(5)
          end
          pause(0.1);
        end

        [obj.y_prev, obj.y_now, obj.y_next] = bc(obj);
        obj.y_next = obj.step();
        for ia=1:length(obj.att_scheme)
          if length(obj.delta >= ia)
            if obj.att_scheme(ia) == 1
              obj.y_next = obj.att1(obj.delta(ia));
            elseif obj.att_scheme(ia) == 2
              obj.y_next = obj.att2(obj.delta(ia));
            elseif obj.att_scheme(ia) == 3
              obj.y_next = obj.att3(obj.delta(ia));
            end
          end
        end

        for mic_num_i = 1:obj.mic_num
          obj.mic(mic_num_i,mic_i) = obj.y_next(obj.mic_pos_i(mic_num_i));
        end

        if any(obj.n == obj.snaps_n)
          obj.snaps(snaps_i,:) = obj.y_next;
          snaps_i = snaps_i + 1;
        end

        obj.n = obj.n + 1;

        obj.y_prev = obj.y_now;
        obj.y_now  = obj.y_next;
        mic_i = mic_i + 1;

        if obj.first
          obj.first = 0;
        end
      end

      if obj.video_file
        close(v);
      end

      for mic_num_i = 1:obj.mic_num
        figure(10+mic_num_i)
        plot(obj.t_axis,obj.mic(mic_num_i,:))
        if obj.mic_plot_png
          print([obj.output_dir, ...
                 obj.output_file_prefix, ...
                 'mic_', num2str(mic_num_i), ...
                 '.png'], '-dpng')
        end
      end
      csvwrite([obj.output_dir, ...
                obj.output_file_prefix, ...
                obj.mic_csvfile], obj.mic');
      if snaps_i >= 2
        figure(40)
        hold on;
        for i = 1:snaps_i -1
          plot(obj.x, obj.snaps(i,:))
          title('snaps')
        end
        hold off
        csvwrite([obj.output_dir, ...
                  obj.output_file_prefix, ...
                  'snaps.csv'],obj.snaps');
        print([obj.output_dir, ...
               obj.output_file_prefix, ...
               'snaps.png'], '-dpng')
      end
    end

    function disp_time(obj, t)
      persistent counter
      if obj.first
        counter = 0;
      end
      counter = counter + 1;
      fprintf('%8.4g',t);
      if counter ~= 10
        fprintf(' | ');
      else
        fprintf('\n');
        counter = 0;
      end
    end

    function y_next = step(obj)

      y_next = obj.y_next;

      for i = [2:length(obj.y_next)-1]
        y_next(i) = (- obj.y_prev(i) ...
                     + 2 * obj.y_now(i) ...
                     + obj.cdtdx2 * (obj.y_now(i-1) ...
                                     - 2 * obj.y_now(i) ...
                                     + obj.y_now(i+1)) );
      end
    end

    function y_next = att1(obj, delta)

      y_next = obj.y_next;

      %% To calculate a third derivarive here D^3_t, it whould have
      %% been necessary to have another time level, say,
      %% y_prev2. Since we dont't have that, we calculate D_t D^2_x,
      %% which is not an additional approximation, but an intermediate
      %% step in obtaining the attenueation D^3_t form first
      %% principles. However, this is not very stable, use of bigger
      %% Nx, or bigger delta breaks the excecution. Something better
      %% is needed.

      for i = [2:length(obj.y_next)-1]
        y_next(i) = ( y_next(i) ...
                      + (0.003*delta*obj.dt/obj.dx^2) ...
                        * (obj.y_now(i+1) ...
                           - 2.0*obj.y_now(i) ...
                           + obj.y_now(i-1) ...
                           - obj.y_prev(i+1) ...
                           + 2.0*obj.y_prev(i) ...
                           - obj.y_prev(i-1)) );
      end
    end

    function y_next = att2(obj, delta)

      y_next = obj.y_next;

      for i = [2:length(obj.y_next)-1]
        y_next(i) = ( y_next(i) ...
                      + (delta*obj.dt) ...
                        * (obj.y_prev(i) - obj.y_now(i)) );
      end
    end

    function y_next = att3(obj, delta)

      y_next = obj.y_next;

      %% Doing nothing
      disp('WARNING: attenuation scheme 3 is not implemented.')
    end

    function D = D_get(obj)
      D = (obj.y_now - obj.y_prev)./obj.dt;
    end

    function [y_left, y_right] = get_split(obj)

      y_left = 0.5*obj.y_next;
      y_right = 0.5*obj.y_next;
      integral = 0;

      for i = [1:length(obj.x)]
        integral = integral + obj.dx * obj.D(i);
        y_left(i) = y_left(i) ...
                    + (1.0/(2.0*obj.c)) * integral;
        y_right(i) = y_right(i) ...
                     - (1.0/(2.0*obj.c)) * integral;
      end
      mean_left = mean(y_left);
      mean_right = mean(y_right);
      mm = 0.5 * (mean_left - mean_right);

      %% This operation ensures mean(y_left) == mean(y_right)
      y_left = y_left - mm;
      y_right = y_right + mm;
    end

    function pause(obj, src, event)

      persistent pause
      if isempty(pause)
        pause = true;
      end
      disp(['pause = ' num2str(pause)])

      if pause
        disp('pause')
        pause = not(pause);
        k = waitforbuttonpress();
        disp(k)
        pause = not(pause);
        disp('release')
      else
        disp('continue')
      end
    end

    function obj = ani_init(obj, figNum, ylim, al_i, title)
      if ishandle(figNum)
        close(figNum)
      end
      figure(figNum)
      if ~obj.is_octave & obj.matlabVerFlag1
        for i=1:length(al_i)
          if i<=length(obj.al_p) & iscell(obj.al_p{i})
            obj.al{al_i(i)} = animatedline(obj.al_p{i}{:});
          else
            obj.al{al_i(i)} = animatedline;
          end
        end
        %% This is not working in octave, DKW
        set(figNum, ...
            'WindowButtonDownFcn', ...
            @(src, event) obj.pause(src, event))
      end
      if obj.is_octave
        set(figNum, ...
            'WindowButtonDownFcn', ...
            @() pause(0.5))
      end
      set(gca, ...
          'xlim', [obj.x(1) obj.x(end)], ...
          'ylim', [ylim(1) ylim(2)]);
      if ~isempty(title)
        if obj.is_octave
          set(gca, ...
              'title', title);
        else
          ax = gca;
          ax.Title.String = title;
        end
      end
      grid on
    end

    function ani_step(obj, figNum, y, ylim, al_i, title)
      if ~obj.is_octave & obj.matlabVerFlag1
        for i=1:length(al_i)
          clearpoints(obj.al{al_i(i)});
          addpoints(obj.al{al_i(i)}, ...
                    obj.x, ...
                    y(:,i));
        end
      else
        figure(figNum);
        for i=1:size(y,2)
            if i==2
              hold on;
            end
            if i<=length(obj.al_p) && iscell(obj.al_p{i})
              plot(obj.x,y(:,i),obj.al_p{i}{:});
            else
              plot(obj.x,y(:,i));
            end
        end
        if ishold
          hold off;
        end
        set(gca, ...
          'xlim', [obj.x(1) obj.x(end)], ...
          'ylim', [ylim(1) ylim(2)]);
        if ~isempty(title)
          set(gca, ...
              'title', title);
        end
        grid on;
        drawnow;
      end
    end
  end
end

%% Local Variables:
%% indent-tabs-mode: nil
%% End:
