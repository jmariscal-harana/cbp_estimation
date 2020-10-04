function [ERROR_POP] = ParameterError(REF,EST,Parameter)

%	Population absolute errors
eval(['[ERROR_POP.Abs.',Parameter,'.Mean, ERROR_POP.Abs.',Parameter,'.SD, ERROR_POP.Abs.',Parameter,'.Max, ERROR_POP.Abs.',Parameter,'.Min] = DatasetError([REF.',Parameter,'],[EST.',Parameter,'],''Abs'');'])

%	Population relative errors
eval(['[ERROR_POP.Rel.',Parameter,'.Mean, ERROR_POP.Rel.',Parameter,'.SD, ERROR_POP.Rel.',Parameter,'.Max, ERROR_POP.Rel.',Parameter,'.Min] = DatasetError([REF.',Parameter,'],[EST.',Parameter,'],''Rel'');'])

end

function [Mean,SD,Max,Min] = DatasetError(REF,EST,ERR_type)
%%	Calculate individual (intrasubject) errors
switch ERR_type
	case 'Abs'	%	Absolute errors
		ERROR		= EST - REF;	%my method
% 		ERROR		= abs(EST - REF);	%jordi's method
	case 'Rel'	%	Relative errors
		ERROR_abs	= EST - REF;	%my method
% 		ERROR_abs	= abs(EST - REF);	%jordi's method
		ERROR		= 100*ERROR_abs./REF;
end

%%	Calculate mean (SD), max and min errors
Mean	= mean(ERROR);
SD		= std(ERROR);
Max		= max(abs(ERROR));
Min		= min(abs(ERROR));

end