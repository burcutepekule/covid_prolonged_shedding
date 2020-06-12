# Set-up

qsum = function(x) c(`50%`=median(x),quantile(x,c(0.025,0.975)))

logit = function(x) log(x/(1-x))

inv.logit = function(x) exp(x)/(1+exp(x))