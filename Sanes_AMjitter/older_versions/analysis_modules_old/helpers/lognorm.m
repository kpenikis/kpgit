function y = lognorm(mu,sigma,x)
y = erfc( (log(x) - mu) / sigma );
end
