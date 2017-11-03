function [] = plotsigmoids(X,pbest,ps,feats)

	D = size(X,2);

	for i=1:D
		subplot(3,3,i);
		r = range(X(:,i));
		xs = linspace( min(X(:,i))-0.05*r, max(X(:,i))+0.05*r, 200)';

		plot(xs, sigmoid(xs, pbest(1+(i-1)*2 : i*2)), 'r-');

		if ~isempty(ps)

			hold on;
			fs = zeros( size(X,1), length(ps));
			for j=1:length(ps)
				plot(xs, sigmoid(xs, ps(j, 1+(i-1)*2 : i*2)), 'k-','color',0.8*[1 1 1]);
				plot(X(:,i), sigmoid(X(:,i), ps(j, 1+(i-1)*2 : i*2)), 'k.','color',0.65*[1 1 1]);
				fs(:,j) = sigmoid(X(:,i), ps(j, 1+(i-1)*2 : i*2));
			end
			plot(xs, sigmoid(xs, pbest(1+(i-1)*2 : i*2)), 'r-');
			hold off;

			refline(1);
			ry = max(fs(:))-min(fs(:));
			ylim([min(min(fs(:))-0.05*ry,min(xs)) max(xs)]);
		end

		xlim([min(xs) max(xs)]);

		title({ strrep(feats{i},'_','\_'), sprintf(['$y=-e^{%.3f} + \\frac{2e^{%.3f}}{1+e^{-xe^{%.3f}}}$'], pbest(i*2), pbest(i*2), pbest(1+(i-1)*2) )},'interpreter','latex');
%		title(feats{i}, 'interpreter','none');
		xlabel('Rosetta score');
		ylabel('Weighted score');
	end
end
