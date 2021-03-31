function map = rwbcm(varargin)
% Red-White-Blue Colormap.
% Example: figure;membrane;colormap(rwbcm(20,1,1,1,'Darken',.9));colorbar;

    p = inputParser;
    validNumber = @(x) mod(x,1) == 0 && isscalar(x);
    p.addOptional('m', 150, validNumber);
    p.addOptional('wshift', 0, validNumber);
    p.addOptional('wspan', 0, validNumber);
    p.addOptional('rbspan', 0, validNumber);
    p.addParameter('Darken', 1, @isscalar);
    p.parse(varargin{:});
    s = p.Results;

	if nargin < 1
	% Set m to the size of the current or default colormap
	% when no argument is given to this function.
	% This is a common preamable in colormap definitions.
	   f = get(groot,'CurrentFigure');
	   if isempty(f)
	      m = size(get(groot,'DefaultFigureColormap'),1);
	   else
	      m = size(f.Colormap,1);
	   end
	end

    % Shift white levels. Tuning: [More Blue](s.m-wspan)/2 <<--- wshift --->> (-s.m+wspan)/2[More Red]
    if s.wshift > (s.m - s.wspan)/2
    	s.wshift = (s.m - s.wspan)/2;
    	warning("That's too Red!");
    elseif s.wshift < (s.wspan - s.m)/2
    	s.wshift = (s.m - s.wspan)/2;
    	warning("That's too Blue!");
    end

    % Pure white levels. Tuning: [Red-Nil White-Blue]0 <<--- s.wspan -->> s.m[Total White]
    if s.wspan < 0
    	s.wspan = 0;
    	warning('White has been extinct!');
    elseif s.wspan > s.m
    	s.wspan = s.m;
    	warning('Colours has been extinct!');
    end

    % Pure red/blue levels. Tuning: [Red-Fades to White-Fades to Blue]0 <<--- s.rbspan --->> Inf[Pure Red-Abrupt change into Pure Blue]
    if s.rbspan < 0
    	s.rbspan = 0;
    	warning('Red-Blue-span must be positive!');
    elseif s.rbspan > s.m/2
    	warning('Fading has been extinct!');
    end

    % The darkening coefficient.
    if s.Darken > 1
    	s.Darken = 1;
    	warning('The darkening coefficient shoult be between 0~1.');
    elseif s.Darken < 0
    	s.Darken = 0;
    	warning('The darkening coefficient shoult be between 0~1.');
    end

    s.m = s.m - s.rbspan*2;

    % Proper definition:
    map = [[linspace(0,1,(s.m-s.wspan)/2-s.wshift),linspace(1,1,s.wspan),linspace(1,1,(s.m-s.wspan)/2+s.wshift)];... % Red Curve
        [linspace(0,1,(s.m-s.wspan)/2-s.wshift),linspace(1,1,s.wspan),linspace(1,0,(s.m-s.wspan)/2+s.wshift)];... % Green Curve
        [linspace(1,1,(s.m-s.wspan)/2-s.wshift),linspace(1,1,s.wspan),linspace(1,0,(s.m-s.wspan)/2+s.wshift),]]'; % Blue Curve
    map = [repmat([0,0,1],[s.rbspan,1]); map; repmat([1,0,0],[s.rbspan,1])];

	% Darken the colormap.
    map = s.Darken*map;

    return;
end
