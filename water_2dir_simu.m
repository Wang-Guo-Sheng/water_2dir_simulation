%% Water 2DIR Simulation (With Cumulant Truncation)
% Example: obj = water_2dir_simu('force.xvg','coord.xvg',figure(),'t2s',0:.1:1.5)

classdef water_2dir_simu
    %% WATER_2DIR_SIMU Generate 2DIR spectrum from calculated frequency trace.
    %   Generate 2DIR spectrum from calculated frequency trace.
    
    properties
    %%
        constants = struct();

        t = [];
        w = [];
        dw = [];
        meanw = NaN;

        tcorr = [];
        FFCF = [];
        D2 = NaN;
        StartPoint = [];
        foo = cfit();
        gof = NaN;
        cv = [];

        g = @(x) x;
        tt = [];
        T1 = [];
        T3 = [];
        R_nr = [];
        R_r = [];
        R_nr_f = [];
        R_r_f = [];
        R = [];

        a = NaN;
        dw1 = NaN;
        w1 = [];

        t2s = [];
        CLS = [];
        StartPoint_CLS = [];
        foo_CLS = cfit();
        cv_CLS = [];
    end

    properties (Constant)
    %%
        fitfunc = @(p, t) (p(1)*cos(p(2)*t).*exp(-p(3)*t) + p(4)*exp(-p(5)*t) + p(6).*exp(-p(7)*t));
        ft = fittype('a1.*cos(k11.*t).*exp(-k12.*t)+a2.*exp(-k2.*t)+a3.*exp(-k3.*t)','independent','t');
        cvperm = [1 4 5 2 6 3 7];

        proto_g = @(p, D2, t) D2*exp(-t.*(p(3)+p(5)+p(7)))/(p(5)^2*p(7)^2*(p(2)^2+p(3)^2)^2) ...
            .*(p(1)*p(5)^2*p(7)^2.*exp((p(5)+p(7)).*t).*(cos(p(2).*t).*(p(3)^2-p(2)^2)-2*p(2)*p(3)*sin(p(2)*t) ...
            + exp(p(3).*t).*(p(2)^2*(p(3).*t+1)+p(3)^2*(p(3).*t-1))) ...
            + (p(2)^2+p(3)^2)^2.*exp(p(3).*t).*(p(6)*p(5)^2*exp(p(5).*t).*(exp(p(7).*t).*(p(7).*t-1)+1) ...
            +p(4)*p(7)^2*exp(p(7).*t).*(exp(p(5).*t).*(p(5)*t-1)+1)));
        nt = 512;
        dtt = 0.0025;
    end

    methods
    %%
        function obj = water_2dir_simu(varargin)
        %% WATER_SIMU 构造此类的实例
            %   此处显示详细说明
            validnum = @(x) isscalar(x) && isnumeric(x);
            validint = @(x) validnum(x) && isinteger(x);
            validsp = @(v) isvector(v) && (length(v) == 7) && isnumeric(v);
            validvec = @(v) isvector(v) && isnumeric(v);

            p = inputParser();
            p.addRequired("force_file_name", @isfile);
            p.addRequired("coord_file_name", @isfile);
            p.addOptional("hf", gcf, @ishandle);
            p.addParameter("nsteps", 100001, validint);
            p.addParameter("t2", 0, validnum);
            p.addParameter("t2s", [], validvec);
            p.addParameter("dt", 0.01, validnum);
            p.addParameter("ncorr", 150, validint);
            p.addParameter("StartPoint", [0.3232, 0.3378, 0.3455, 30, 17, 8, 1.9], validsp);
            p.parse(varargin{:});
            s = p.Results;

            obj.constants = obj.get_constants();
            obj.t = (0:s.nsteps - 1)*s.dt;
            obj.w = NaN*ones(size(obj.t));
            obj.dw = obj.w;
            obj.tcorr = (0:s.ncorr - 1)*s.dt;
            if isempty(s.t2s)
                s.t2s = obj.tcorr;
            end
            obj.FFCF = NaN*ones(size(obj.tcorr));
            obj.StartPoint = s.StartPoint;
            obj.tt = (0:obj.nt - 1)*obj.dtt;
            [obj.T1, obj.T3] = meshgrid(obj.tt, obj.tt);
            obj.t2s = s.t2s;

            obj = obj.read_freq(s.force_file_name, s.coord_file_name);
            obj = obj.compute_FFCF(s.ncorr);

            ax1 = subplot(121, 'Parent', s.hf);
            obj.plot_FFCF(ax1);
            obj = obj.fit_FFCF(ax1);


            obj = obj.compute_spectrum(s.t2);

            ax2 = subplot(122, 'Parent', s.hf);
            obj.plot_spectrum(ax2);
            title(ax2, '2D LineShape');

            obj.CLS = obj.compute_CLS(obj);

            obj.plot_CLS(ax1);

            % obj.StartPoint_CLS = obj.StartPoint*obj.CLS(1);
            % obj = obj.fit_CLS(ax1);

            title(ax1, 'FFCF and CLS');
            legend(ax1, 'show');

            return;
        end
        

        function obj = read_freq(obj, force_file_name, coord_file_name)
        %%
           [force_O, force_H1] = obj.read_xvg(force_file_name);
           [pos_O, pos_H1] = obj.read_xvg(coord_file_name);
           [obj.w, obj.dw, obj.meanw] = obj.get_freq(force_O, force_H1, pos_O, pos_H1);
           return;
        end


        function hf = plot_trajectory(obj, ax)
        %%
            plot(ax, obj.t, obj.w);
            xlabel(ax, 't/ps');
            ylabel(ax, '\omega/cm^{-1}');
        end


        function obj = compute_FFCF(obj, ncorr)
        %%
            for idx = 0:ncorr - 1
                corr_itr = obj.dw.*circshift(obj.dw, -idx);
                obj.FFCF(idx + 1) = mean(corr_itr(1:length(obj.dw) - idx));
            end
            obj.D2 = obj.FFCF(1);
            obj.FFCF = obj.FFCF/obj.D2;

            return;
        end


        function hf = plot_FFCF(obj, ax)
        %%
            scatter(ax, obj.tcorr, obj.FFCF, 'k.', 'DisplayName', 'FFCF data');
            xlabel(ax, 't/ps');
            ylabel(ax, 'Normalized FFCF');
        end


        function [obj, hf] = fit_FFCF(obj, ax)
        %%
            ih = ishold(ax);
            if ~ih
                hold(ax, 'on');
            end

            x = obj.tcorr';
            y = obj.FFCF';

            plot(ax, x, obj.fitfunc(obj.StartPoint(obj.cvperm), x), 'b-', 'DisplayName', 'Start Point');

            obj.foo = fit(x, y, obj.ft, ...
                'StartPoint', obj.StartPoint, ...
                'Lower', [0, 0, 0, 1e1, 1e0, 1e0, 1e-1], ...
                'Upper', [1, 1, 1, 1e2, 1e3, 1e2, 1e1]);
            obj.cv = coeffvalues(obj.foo);
            obj.gof = goodnessOfFit(obj.foo(x), y, 'MSE');

            %%
            MAXIT = 100;
            for idx = 1:MAXIT
                obj.foo = fit(x, y, obj.ft, ...
                    'StartPoint', obj.cv, ...
                    'Lower', [0, 0, 0, 1e1, 1e0, 1e0, 1e-1], ...
                    'Upper', [1, 1, 1, 1e2, 1e3, 1e2, 1e1]);
                obj.cv = coeffvalues(obj.foo);

                gof = goodnessOfFit(obj.foo(x), y, 'MSE');
                if obj.gof == gof
                    break;
                else
                    obj.gof = gof;
                end
            end

            %%
            plot(ax, x, obj.foo(x), 'r-', 'DisplayName', 'Fit');

            obj.g = @(t) obj.proto_g(obj.cv(obj.cvperm), obj.D2, t);

            if ~ih
                hold(ax, 'off');
            end
            return;
        end


        function obj = compute_spectrum(obj, t2)
        %%
            const = obj.get_constants();

            [R, obj] = obj.compute_spectrum_static(obj, t2);
            obj.R = R;
    
            obj.a = 1/obj.dtt/(const.c*1e-12);
            obj.dw1 = obj.a/(2*length(obj.tt));
            obj.w1 = (-obj.a/2+obj.dw1:obj.dw1:obj.a/2);
        end


        function plot_spectrum(obj, varargin)
        %%
            validnum = @(x) isscalar(x) && isnumeric(x);
            validrg = @(v) isvector(v) && (length(v) == 2) && isnumeric(v);

            p = inputParser();
            p.addRequired("ax", @ishandle);
            p.addParameter("n_contours", 20, validnum);
            p.addParameter("plot_range", [2800 3600], validrg);
            p.parse(varargin{:});
            s = p.Results;
    
            ind1 = find(obj.w1>=s.plot_range(1) & obj.w1<=s.plot_range(2));
            contourf(s.ax, obj.w1(ind1), obj.w1(ind1), obj.R(ind1,ind1)',s.n_contours);

            axis(s.ax,'equal','tight')
            % set(s.ax,'xtick',[s.plot_range(1):1000:s.plot_range(end)])
            % set(s.ax,'ytick',[s.plot_range(1)-1:1000:s.plot_range(end)])
            line(s.ax,[s.plot_range(1) s.plot_range(end)], [s.plot_range(1) s.plot_range(end)],'Color','k');

            xlabel(s.ax,'\omega_3/cm^{-1}');
            ylabel(s.ax,'\omega_1/cm^{-1}');

            try
                colormap(s.ax, 'rwbcm');
            catch
            end
        end


        function obj = plot_CLS(obj, ax)
        %%
            ih = ishold(ax);
            if ~ih
                hold(ax, 'on');
            end
            clss = plot(ax, obj.t2s, obj.CLS, 'm.', 'DisplayName', 'CLS data');
            xlabel(ax, 't/ps');
            ylabel(ax, 'CLS');
            if ~ih
                hold(ax, 'off');
            end
        end


        function [obj, hf] = fit_CLS(obj, ax)
        %%
            ih = ishold(ax);
            if ~ih
                hold(ax, 'on');
            end
            plot(ax, obj.t2s, obj.fitfunc(obj.StartPoint_CLS(obj.cvperm), obj.t2s), 'b-', 'DisplayName', 'Start Point');

            obj.foo_CLS = fit(obj.t2s', obj.CLS', obj.ft, ...
                'StartPoint', obj.StartPoint_CLS, ...
                'Lower', [0, 0, 0, 1e1, 1e0, 1e0, 1e-1], ...
                'Upper', [1, 1, 1, 1e2, 1e3, 1e2, 1e1]);
            plot(ax, obj.t2s, obj.foo(obj.t2s), 'r-', 'DisplayName', 'Fit');

            obj.cv_CLS = coeffvalues(obj.foo);

            if ~ih
                hold(ax, 'off');
            end
            return;
        end

        function outputArg = method1(obj,inputArg)
            %METHOD1 此处显示有关此方法的摘要
            %   此处显示详细说明
            outputArg = obj.Property1 + inputArg;
        end
    end

    methods (Static)
        function constants = get_constants()
        %%
            constants = struct();
            constants.N_A = 6.02214086e23;
            constants.c = 2.99792458e10;    % cm/s
            constants.h = 6.62607004e-34;
            constants.mass = 1.66053906660e-27;
            constants.Delta = 240; % Anharmonic Shift.
            constants.Delta_ps = constants.Delta*(2*pi*constants.c*1e-12);
            constants.w_0 = 3550;
            constants.mH = 1.00794;
            constants.mO = 15.999;
            constants.red_mass = constants.mH*constants.mO/(constants.mH + constants.mO);
            constants.a = -1.5*sqrt(constants.Delta/constants.c/constants.h/constants.red_mass/constants.mass)/(2*pi*constants.c*1e-12)/constants.w_0/constants.N_A;
        end

        function [O_data, H1_data] = read_xvg(fname)
        %%
            fid =fopen(fname);
            data = fscanf(fid,'%f\t');
            fclose(fid);

            data = reshape(data,[7,length(data)/7]).';
            O_data = data(:,2:4);
            H1_data = data(:,5:7);

            return;
        end

        function varargout = get_freq(varargin)
        %%
            nargoutchk(0, 3);
            p = inputParser();
            p.addRequired('force_O', @isnumeric);
            p.addRequired('force_H1', @isnumeric);
            p.addRequired('pos_O', @isnumeric);
            p.addRequired('pos_H1', @isnumeric);
            p.parse(varargin{:});
            s = p.Results();

            %%
            const = water_2dir_simu.get_constants();

            bond = s.pos_H1 - s.pos_O;
            bond_length = sqrt(sum(bond.^2,2));
            u_bond = bond./bond_length;
            force = (s.force_H1/const.mH - s.force_O/const.mO)*const.red_mass;
            force_proj = sum(force.*u_bond, 2);

            w = const.a*force_proj + const.w_0;
            meanw = mean(w);
            dw = w - meanw;

            dw = dw*(2*pi*const.c*1e-12);
            meanw = meanw*(2*pi*const.c*1e-12);

            %%
            varargout = {w, dw, meanw};
            return;
        end


        function [R, obj] = compute_spectrum_static(obj, t2)
        %%
            const = obj.get_constants();

            obj.R_nr = exp(-obj.g(obj.T1)-obj.g(t2)-obj.g(obj.T3)+obj.g(obj.T1+t2)+obj.g(t2+obj.T3)-obj.g(obj.T1+t2+obj.T3))...
            .*(exp(-1i.*(obj.T3+obj.T1)*obj.meanw)-exp(-1i.*(obj.T3.*(obj.meanw-const.Delta_ps)+obj.T1*obj.meanw)));
            obj.R_r = exp(-obj.g(obj.T1)+obj.g(t2)-obj.g(obj.T3)-obj.g(obj.T1+t2)-obj.g(t2+obj.T3)+obj.g(obj.T1+t2+obj.T3))...
            .*(exp(-1i.*(obj.T3-obj.T1)*obj.meanw)-exp(-1i.*(obj.T3.*(obj.meanw-const.Delta_ps)-obj.T1*obj.meanw)));

            obj.R_r(:,1) = obj.R_r(:,1)./2;
            obj.R_r(1,:) = obj.R_r(1,:)./2;
            obj.R_nr(:,1) = obj.R_nr(:,1)./2;
            obj.R_nr(1,:) = obj.R_nr(1,:)./2;

            obj.R_nr_f = fftshift(ifft2(obj.R_nr,2*obj.nt,2*obj.nt));
            obj.R_r_f = fftshift(ifft2(obj.R_r, 2*obj.nt,2*obj.nt));
            obj.R_r_f = fliplr(circshift(obj.R_r_f,[0 -1]));
            obj.R = real(obj.R_r_f + obj.R_nr_f);

            obj.R_nr_f = fftshift(ifft2(obj.R_nr,2*obj.nt,2*obj.nt));
            obj.R_r_f = fftshift(ifft2(obj.R_r, 2*obj.nt,2*obj.nt));
            obj.R_r_f = fliplr(circshift(obj.R_r_f,[0 -1]));
            R = real(obj.R_r_f + obj.R_nr_f);
    
            return;
        end


        function cls = max_line_slope(R, ind1, hw)
        %%
            R = R(ind1, ind1);
            hw = fix(hw);
            
            mlx = 1:length(R);
            [~,mly] = water_2dir_simu.symmin(-R);
            mlx = mlx(2:end-1);
            mly = mly(2:end-1);
            
            [~, I] = max(R(:));
            [r, c] = ind2sub(size(R),I);

            idx = -hw:hw;

            p = polyfit(mlx(r + idx), mly(c + idx)', 1);
            cls = p(1);
        end


        function [M, I] = symmin(v)
        %% Get minimum value and position, while returning a center index if the
        % minimum is a flat plateau.
            sz = size(v);
            if isvector(v)
                [~,d] = max(sz);
            elseif ismatrix(v)
                d = 2; % Find by row, i.e. along the probe axis.
            else
                error('`symmovmean` is not intended to handle arrays with dimensions other than 1 or 2.');
            end

            [M,I] = min(v,[],d);
            [~,flI] = min(flip(v,d),[],d);
            flI = sz(d) - flI + 1;
            I = round(.5*(I+flip(flI,d)));
            return;
        end


        function clss = compute_CLS(obj, varargin)
        %%
            validnum = @(x) isscalar(x) && isnumeric(x);
            validrg = @(v) isvector(v) && (length(v) == 2) && isnumeric(v);
            p = inputParser();
            p.addParameter("peak_range", [2800 3600], validrg);
            p.addParameter("half_width", 100, validnum);
            p.parse(varargin{:});
            s = p.Results;

            ind1 = (obj.w1>=s.peak_range(1) & obj.w1<=s.peak_range(2));

            func = @(t2) obj.max_line_slope(obj.compute_spectrum_static(obj, t2), ind1, s.half_width/obj.dw1);
            clss = arrayfun(func, obj.t2s);
            return;
        end
    end
end

