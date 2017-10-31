
% sigmoidal function with equal upper and lower bounds (gives f(0)=0)
function f = sigmoid(X,p)

	a = p(1:2:end);
	r = p(1+1:2:end);

	f = sum( -exp(r) + 2*exp(r) ./ (1+exp(-exp(a).*X)),2);

%	logitlin = @(X, a,r) sum( -exp(r) + 2*exp(r) ./ (1+exp(-exp(a).*X)),2);
%	logitlinp = @(X,p) logitlin(X, p(1:2:end), p(1+1:2:end));
	
end

