% load ../data/zemu.mat
%
% fref = zbrr.ros;
% ftal = zbrt.ros;
% fcon = zc.ros;
%
% %% backrub-ref
%
% Rref = sigmoidfit(zbrr.X,zbrr.y,50,1000);
%
% plotsigmoids(zbrr.X,Rref.phat, Rref.ps(1:100,:), zbrr.feats);
% print('zemu_sigmoid2_ref_feats.png','-dpng','-r300');
%
% plotsample(Rref);
% print('zemu_sigmoid2_ref_posterior.png','-dpng','-r300');


%% backrub-talaris

talaris_table = readtable('zemu-backrub-1.2-50-30000-t14.csv');
talaris_fields = {'fa_sol', 'hbond_sc', 'hbond_bb_sc', 'fa_rep', 'fa_elec', 'hbond_lr_bb', 'fa_atr'};
[m,n] = size(talaris_fields)
pred_data = zeros( 1240, n );
exp_data = zeros( 1240, 1 );
for i = 1:n
    field_name = char(talaris_fields(i));
    pred_data(:,i) = talaris_table.(field_name);
    %% exp_data(:,i) = talaris_table.ExperimentalDDG(i);
end
exp_data = talaris_table.ExperimentalDDG;

Rtal = sigmoidfit(pred_data, exp_data, 50, 1000);

plotsigmoids(pred_data, Rtal.phat, Rtal.ps(1:100,:), talaris_fields);
print('zemu_sigmoid2_tal_feats.png','-dpng','-r300');

plotsample(Rtal);
print('zemu_sigmoid2_tal_posterior.png','-dpng','-r300');

return

%% control

Rcon = sigmoidfit(zc.X,zc.y,50,1000);

plotsigmoids(zc.X,Rcon.phat, Rcon.ps(1:100,:), zc.feats);
print('zemu_sigmoid2_con_feats.png','-dpng','-r300');

plotsample(Rcon);
print('zemu_sigmoid2_con_posterior.png','-dpng','-r300');

%% plot correlations from all

subplot(321);
plot(zc.ros, zc.y,'r.');
title(sprintf('Control, corr %.3f, MAE %.3f, MSE %.3f', corr(zc.ros,zc.y), mean(abs(zc.ros-zc.y)),mean((zc.ros-zc.y).^2) ));
refline(1); xlim([-4 10]); ylim([-7 10]);


subplot(323);
plot(zbrt.ros, zbrt.y,'r.');
title(sprintf('Backrub Talaris, corr %.3f, MAE %.3f, MSE %.3f', corr(zbrt.ros,zbrt.y), mean(abs(zbrt.ros-zbrt.y)),mean((zbrt.ros-zbrt.y).^2) ));
refline(1); xlim([-4 10]); ylim([-7 10]);


subplot(325);
plot(zbrr.ros, zbrr.y,'r.');
title(sprintf('Backrub REF, corr %.3f, MAE %.3f, MSE %.3f', corr(zbrr.ros,zbrr.y), mean(abs(zbrr.ros-zbrr.y)),mean((zbrr.ros-zbrr.y).^2) ));
refline(1); xlim([-4 10]); ylim([-7 10]);


subplot(322);
plot(Rcon.fhat, zc.y,'r.');
hold on;
plot( [min(Rcon.fs')' max(Rcon.fs')']', [zc.y zc.y]', '-', 'color', 0.65*[1 1 1]);
plot(Rcon.fhat, zc.y,'r.');
hold off;

title(sprintf('Control, corr %.3f, MAE %.3f, MSE %.3f', corr(Rcon.fhat,zc.y), mean(abs(Rcon.fhat-zc.y)),mean((Rcon.fhat-zc.y).^2) ));
refline(1); xlim([-4 10]); ylim([-7 10]);


subplot(324);
plot(Rtal.fhat, zbrt.y,'r.');
hold on;
plot( [min(Rtal.fs')' max(Rtal.fs')']', [zbrt.y zbrt.y]', '-', 'color', 0.65*[1 1 1]);
plot(Rtal.fhat, zbrt.y,'r.');
hold off;

title(sprintf('Backrub Talaris, corr %.3f, MAE %.3f, MSE %.3f', corr(Rtal.fhat,zbrt.y), mean(abs(Rtal.fhat-zbrt.y)),mean((Rtal.fhat-zbrt.y).^2) ));
refline(1); xlim([-4 10]); ylim([-7 10]);

subplot(326);
plot(Rref.fhat, zbrr.y,'r.');
title(sprintf('Backrub REF, corr %.3f, MAE %.3f, MSE %.3f', corr(Rref.fhat,zbrr.y), mean(abs(Rref.fhat-zbrr.y)),mean((Rref.fhat-zbrr.y).^2) ));
hold on;
plot( [min(Rref.fs')' max(Rref.fs')']', [zbrr.y zbrr.y]', '-', 'color', 0.65*[1 1 1]);
plot(Rref.fhat, zbrr.y,'r.');
hold off;
refline(1); xlim([-4 10]); ylim([-7 10]);

print('zemu_sigmoid2_corrs.png','-dpng','-r300');




%% backrub-talaris 10-CV

folds = cvpartition(1240,'Kfold',10);
fs = zeros(1240,1);
for f=1:10
	tr = folds.training(f);
	ts = folds.test(f);

	% train
	Rf = sigmoidfit(zbrt.X(tr,:),zbrt.y(tr));
	% test
	fs(ts) = sigmoid(zbrt.X(ts,:), Rf.phat);
end


Rtal = sigmoidfit(zbrt.X,zbrt.y,50,1000);
plotsigmoids(zbrt.X,Rtal.phat, Rtal.ps(1:100,:), zbrt.feats);
plotsample(Rtal);
