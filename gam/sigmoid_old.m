
% three datasets:
% - zbrr : 8 scores with REF weights      [backrub]
% - zbrt : 7 scores with talaris weights  [backrub]
% - zc   : 7 scores with talaris weights  [no backrub]
load ../data/zemu.mat

X = zbrr.X;
y = zbrr.y;
feats = zbrr.feats;
fref = zbrr.ros;
ftal = zbrt.ros;
fdef = zc.ros;

mae = @(x,y) 


logitlin = @(X, a,b,c) sum(2*c .* (-b + 2*b ./ (1+exp(-a.*X))),2);
logitlinp = @(X,p) logitlin(X, p(1:7), p(8:14), p(15:21));
rmse = @(X,p,y) mean(( logitlinp(X,p)-y).^2);
p0 = [0.1*ones(1,7) 10*ones(1,7) zeros(1,7)];

%% ML solution
options = optimoptions('fminunc','MaxFunctionEvaluations',100000,'MaxIterations',10000);
phat = fminunc( @(p) rmse(X,p,y), p0,options);
fhat = logitlinp(X,phat);

% plot phat
for i=1:7
	subplot(3,3,i);
	xs = linspace( min(X(:,i))-1, max(X(:,i))+1, 200)';
	
	plot(xs, logitlin(xs, phat(i),phat(i+7),phat(i+14)), 'r-');
	hold on;
	plot(X(:,i), logitlin(X(:,i), phat(i),phat(i+7),phat(i+14)), 'k.');
	hold off;
	
	refline(1);
	xlim([min(xs), max(xs)]);
	
	title(feats{i},'interpreter','none');
	xlabel('Rosetta score');
	ylabel('Weighted score');
end



%% plot corr
figure(1);
subplot(131);
plot(zbrt.ros, y,'.');
title(sprintf('Talaris, corr %.3f, MAE %.3f, MSE %.3f', corr(zbrt.ros,y),   ));
refline(1);

subplot(132);
plot(zbrr.ros, y,'.');
title(sprintf('REF, corr %.3f', corr(zbrr.ros,y)));
refline(1);

subplot(133);
plot(fhat,y,'.');
title(sprintf('Sigmoid, corr %.3f', corr(fhat,y)));
refline(1);



%% posterior sampling
delta = 2;

proppdf = @(p1,p2) prod(unifpdf(p2-p1,-delta,delta));
proprnd = @(p) p + 0.001*randn(1,21);

logpdf = @(p) logmvnpdf(logitlinp(X,p),y,0.5);

M = 50000;
%psample = mhsample(p0, M, 'logpdf',logpdf, 'proprnd',proprnd,'proppdf',proppdf,'thin',5,'burnin',1000);
psample = mhsample(phat, M, 'logpdf',logpdf, 'proprnd',proprnd,'proppdf',proppdf,'thin',5,'burnin',0);

ps = zeros(M,1);
mae = zeros(M,1);
mse = zeros(M,1);
cc = zeros(M,1);
for i=1:M
	ps(i) = logpdf(psample(i,:));
	mae(i) = mean(abs( logitlinp(X,psample(i,:))-y));
	mse(i) = mean(( logitlinp(X,psample(i,:))-y).^2);
	cc(i) = corr(logitlinp(X,psample(i,:)),y);
end

%% plot samples
figure(2);
subplot(411); plot(ps); title('logpdf');
subplot(412); plot(cc); title('Corr');
subplot(413); plot(mse); title('MSE');
subplot(414); plot(mae); title('MAE');

%print('zemu_sigmoid_solutions.png','-dpng','-r300');

[v,I] = sort(mae,'ascend');
I = randperm(M);
pmode = psample(I(1),:);

%% plot posterior
for i=1:7
	subplot(3,3,i);
	xs = linspace( min(X(:,i))-1, max(X(:,i))+1, 200)';
	
	plot(xs, logitlin(xs, phat(i),phat(i+7),phat(i+14)), 'r-');
	hold on;
	for j=1:100
		plot(xs, logitlin(xs, psample(I(j),i),psample(I(j),i+7),psample(I(j),i+14)), 'k-','color',0.8*[1 1 1]);
		plot(X(:,i), logitlin(X(:,i), psample(I(j),i),psample(I(j),i+7),psample(I(j),i+14)),'k.','color',0.65*[1 1 1]);
	end
	plot(xs, logitlin(xs, phat(i),phat(i+7),phat(i+14)), 'r-','linewidth',1.5);
	hold off;

	refline(1);
	xlim([min(xs), max(xs)]);
	ylim([min(xs) max(xs)]);
	
	title({feats{i}, sprintf('%.2f / (1 + e^%.2f x) - %.2f', phat(i+14)*2*phat(i+7), phat(i),phat(i+7)*phat(i+14) )},'interpreter','none');
	
	xlabel('Rosetta score');
	ylabel('Weighted score');
end

print('zemu_sigmoid_feats.png','-dpng','-r300');

figure(2);
subplot(131); histogram(cc(I(1:100))); title('Corr');
subplot(132); histogram(mae(I(1:100))); title('MAE');
subplot(133); histogram(mse(I(1:100))); title('MSE');

%% plot corrs
figure(1);
subplot(131);
plot(zbrt.ros, y,'.');
title(sprintf('Talaris, corr %.3f', corr(zbrt.ros,y)));
refline(1);

subplot(132);
plot(zbrr.ros, y,'.');
title(sprintf('REF, corr %.3f', corr(zbrr.ros,y)));
refline(1);

subplot(133);
plot(fhat,y,'.');
title(sprintf('Sigmoid, corr %.3f', corr(fhat,y)));
refline(1);

fs = zeros(1240,100);
for j=1:100
	fs(:,j) = logitlinp(X,psample(I(j),:));
end
errorbar(fhat,y,[],[],sqrt(var(fs'))',sqrt(var(fs'))','.');
refline(1);




% plot top posterior samples
for i=1:6
	subplot(2,3,i);
	xs = linspace( min(X(:,i))-1, max(X(:,i))+1, 200)';
	
%	plot(xs, logitlin(xs, phat(i),phat(i+6),phat(i+12)), 'r-');
	hold on;
	for j=1:4
%		plot(xs, logitlin(xs, psample(I(j),i),psample(I(j),i+6),psample(I(j),i+12)), 'k-','color',0.8*[1 1 1]);
		plot(X(:,i), logitlin(X(:,i), psample(I(j),i),psample(I(j),i+6),psample(I(j),i+12)),'.');
	end
%	plot(xs, logitlin(xs, phat(i),phat(i+6),phat(i+12)), 'r-','linewidth',1.5);
	hold off;

	refline(1);
	xlim([min(xs), max(xs)]);
	ylim([min(xs) max(xs)]);
	
	title({feats{i}, sprintf('%.2f / (1 + e^%.2f x) - %.2f', phat(i+12)*2*phat(i+6), phat(i),phat(i+6)*phat(i+12) )},'interpreter','none');
	
	xlabel('Rosetta score');
	ylabel('Weighted score');
end
print('zemu_sigmoid_feats_ex.png','-dpng','-r300');



% individual terms
ts = zeros(1105,6,100);
for i=1:6
	for j=1:100
		ts(:,i,j) = logitlin(X(:,i), psample(I(j),i),psample(I(j),i+6),psample(I(j),i+12));
	end
end

for i=1:6
	for j=1:6
		subplot(6,6,(i-1)*6+j);
		plot(ts(:,i,2), ts(:,j,2), '.'); 
		title('Post-transform score');
		xlabel(feats{i},'interpreter','none');
		ylabel(feats{j},'interpreter','none');
		xlim([min(ts(:,i,1)) max(ts(:,i,1))]);
		ylim([min(ts(:,i,1)) max(ts(:,i,1))]);
	end
end



% plot correlations
figure(1);
subplot(121);
plin = sum(X,2);
plot(plin, y, 'k.');
title(sprintf('Default, corr %.3f mae %.3f', corr(plin,y), mean(abs(plin-y))));
refline(1);
xlabel('Score');
ylabel('Experimental');

subplot(122);

plog = logitlinp(X,pmode);

plot(plog, y, 'k.');
hold on;
for j=1:100
	plot(logitlinp(X,psample(I(j),:)), y, 'k.','color',0.7*[1 1 1]);
end
plot(plog, y, 'k.');
hold off;
title(sprintf('Logistic, corr %.3f mae %.3f', corr(plog,y), mean(abs(plog-y))));
xlim([-4 10]);
xlabel('Score');
ylabel('Experimental');
refline(1);


print('zemu_sigmoid_corr.png','-dpng','-r300');


%% plot corr-diffs
figure();
plog = logitlinp(X,pmode);

plot(plog, y, 'k.');
hold on;
for j=1:100
	plot(logitlinp(X,psample(I(j),:)), y, 'k.','color',0.7*[1 1 1]);
end
plot(plog, y, 'k.');

Ip = abs(plog-y) < abs(plin-y);
In = abs(plog-y) > abs(plin-y);

plot([plog(Ip) plin(Ip)]', [y(Ip) y(Ip)]', 'r-');
plot([plog(In) plin(In)]', [y(In) y(In)]', 'b-');


%plot([logitlinp(X,pmode) plin]', [y y]', 'b-');
hold off;
title(sprintf('Logistic, corr %.3f mae %.3f', corr(plog,y), mean(abs(plog-y))));
refline(1);
xlabel('Score');
ylabel('Experimental');


print('zemu_sigmoid_corr_diffs.png','-dpng','-r300');







plot(xs, logitlinp(xs,p0), xs, logitlinp(xs,phat));



logit = @(X, a,b,fmax,fmin) fmin + (fmax-fmin) ./ (1+exp(-a.*(X-b)));
func = @(X,p) logit(X, p(1:6), p(7:12), p(13:18), p(19:24));


logit = @(X,a,b,fmax,fmin) sum( fmin + (fmax-fmin) ./ (1+exp(-a.*(X-b))),2);
func = @(X,p) logit(X, p(1:6), p(7:12), p(13:18), p(19:24));

cor_ros = corr(sum(X,2), y);

%% ML solution
p0 = [-0.7*ones(1,6) zeros(1,6) max(X) min(X)];

options = optimoptions('fminunc','MaxFunctionEvaluations',100000,'MaxIterations',10000);
phat = fminunc( @(p) mean(( func(X,p) - y).^2), p0,options);

%% posterior
proppdf = @(p1,p2) prod(unifpdf(p2-p1,-1,1));
proprnd = @(p) p + 0.01*randn(1,24);
%logpdf = @(p) sum(log(normpdf( func(X,p), y, 0.3)));
logpdf = @(p) logmvnpdf(func(X,p),y,0.3);

psample = mhsample(p0, 5000, 'logpdf',logpdf, 'proprnd',proprnd,'proppdf',proppdf,'thin',5,'burnin',5000);

ps = zeros(5000,1);
for i=1:5000
	ps(i) = logpdf(psample(i,:));
end
plot(ps); 


I = find(ps > quantile(ps,0.95));
nnz(I)


%% plot

figure(1);
subplot(121);
plin = sum(X,2);
plot(plin, y, 'k.');
title(sprintf('Default, corr %.3f mse %.3f', corr(plin,y), mean((plin-y).^2)));
refline(1);

subplot(122);
plog = func(X, phat(1:6),phat(7:12),phat(13:18) );
plot(plog, y, 'k.');
title(sprintf('Logistic, corr %.3f mse %.3f', corr(plog,y), mean((plog-y).^2)));
refline(1);


[i,v] = max(ps);
phat = psample(v,:);

figure(2);
for i=1:6
	subplot(2,3,i);
	xs = linspace( min(X(:,i))-1, max(X(:,i))+1, 200)';
	
	plot(xs, logit(xs, phat(i),phat(i+6),phat(i+12),phat(i+18)), 'r-');
	hold on;
	
	for j=1:nnz(I)
%		plot(X(:,i), func(X(:,i), psample(I(j),i),psample(I(j),i+6),psample(I(j),i+12) ), 'k.', 'color',0.6*[1 1 1 ]);		
		plot(X(:,i), logit(X(:,i), psample(I(j),i), psample(I(j),i+6), psample(I(j),i+12),psample(I(j),i+18)), 'k.', 'color',0.6*[1 1 1 ]);		
	end
	plot(xs, logit(xs, phat(i),phat(i+6),phat(i+12),phat(i+18)), 'r-');

	hold off;
	refline(1);
	xlim([min(xs), max(xs)]);
	
	title(feats{i},'interpreter','none');
	xlabel('Rosetta score');
	ylabel('Weighted score');
end

figure(3);
plog = func(X, phat(1:6),phat(7:12),phat(13:18) );
plot(plog,y,'k.');
hold on;
for j=1:nnz(I)
	plot(func(X, psample(I(j),1:6),psample(I(j),7:12),psample(I(j),13:18) ), y, 'k.', 'color',0.6*[1 1 1 ]);
end
plot(plog,y,'k.');
hold off;
refline(1)


xs = linspace( min(X(:,1))-1, max(X(:,1))+1, 200)';
plot(xs, func(xs, phat(1),phat(1+6),phat(1+12)), 'r-');
hold on;

for i=1:175
	plot(X(:,1), func(X(:,1), psample(I(i),1),psample(I(i),1+6),psample(I(i),1+12) ), 'k.', 'color',0.6*[1 1 1 ]);
end
plot(xs, func(xs, phat(1),phat(1+6),phat(1+12)), 'r-', 'linewidth',2);
hold off;

refline(1);
title(feats{i},'interpreter','none');
xlabel('Rosetta score');
ylabel('Weighted score');




