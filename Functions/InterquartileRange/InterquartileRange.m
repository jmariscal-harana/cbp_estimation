function [Q1,Q2,Q3] = InterquartileRange (Data)
Data = sort(Data);

n = length(Data);
if n < 2
	error('One single data point provided: check vector size')
end

if mod(n,2) == 0
	Q1 = median(Data(1:n/2));
	Q3 = median(Data(n/2+1:n));
elseif mod(n,2) == 1
	Q1 = median(Data(1:(n+1)/2));
	Q3 = median(Data((n+1)/2:n));
else
	error('What?!')
end

Q2 = median(Data);

end