function [] = plotsample(res)
	
	subplot(221); histogram(res.lp); title('logpdf');	
	
	subplot(222); histogram(res.cc); title('Corr');
	
	subplot(223); histogram(res.mse); title('MSE');
	
	subplot(224); histogram(res.mae); title('MAE');

end
