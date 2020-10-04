function [Wave1,Wave2,Time_out] = InterpolateSpline(Time_in,Wave1,Wave2)

if length(Wave1) > length(Wave2)
	%upsample Wave2
	t_short		= linspace(Time_in(1),Time_in(end),length(Wave2));
	t_long		= linspace(Time_in(1),Time_in(end),length(Wave1));
	Wave2		= spline(t_short,Wave2,t_long);
	Time_in		= t_long; 
elseif length(Wave1) < length(Wave2)
	%upsample Wave1
	t_short		= linspace(Time_in(1),Time_in(end),length(Wave1));
	t_long		= linspace(Time_in(1),Time_in(end),length(Wave2));
	Wave1		= spline(t_short,Wave1,t_long);
	Time_in		= t_long;
end

Wave1 = Wave1(:);
Wave2 = Wave2(:);
Time_out = Time_in(:);

end