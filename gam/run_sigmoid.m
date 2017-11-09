% load ../data/zemu.mat
%
% fref = zbrr.ros;
% ftal = zbrt.ros;
% fcon = zc.ros;

talaris_table = readtable('../data/by_step/zemu_1.2-60000_rscript_simplified-t14-id_50-30000-partial.csv');
ref_table = readtable('../data/by_step/zemu-backrub-1.2-50-30000-REF-v2.csv');
control_table = readtable('../data/by_step/zemu_control-69aa526-id_50-00008-partial.csv');
talaris_fields = {'fa_sol', 'hbond_sc', 'hbond_bb_sc', 'fa_rep', 'fa_elec', 'hbond_lr_bb', 'fa_atr'};
ref_fields = {'fa_sol', 'hbond_sc', 'hbond_bb_sc', 'fa_rep', 'fa_elec', 'hbond_lr_bb', 'fa_atr', 'lk_ball_wtd'};

%% backrub-ref

[m,n] = size(ref_fields)
ref_pred_data = zeros( 1240, n );
for i = 1:n
    field_name = char(ref_fields(i));
    ref_pred_data(:,i) = ref_table.(field_name);
end
exp_data = ref_table.ExperimentalDDG;

Rref = sigmoidfit(ref_pred_data, exp_data, 50, 1000);

plotsigmoids(ref_pred_data, Rref.phat, Rref.ps(1:100,:), ref_fields);
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 12 10];
print('zemu_sigmoid2_ref_feats.png','-dpng','-r400');

plotsample(Rref);
print('zemu_sigmoid2_ref_posterior.png','-dpng','-r300');

fit_terms = fitterms(ref_pred_data, Rref.phat, ref_fields, exp_data);
writetable( fit_terms, 'ref_GAM_terms.csv' )


%% backrub-talaris

[m,n] = size(talaris_fields)
tal_pred_data = zeros( 1240, n );
for i = 1:n
    field_name = char(talaris_fields(i));
    tal_pred_data(:,i) = talaris_table.(field_name);
end
assert( isequal( exp_data, talaris_table.ExperimentalDDG ) );
exp_data = talaris_table.ExperimentalDDG;

Rtal = sigmoidfit(tal_pred_data, exp_data, 50, 1000);

plotsigmoids(tal_pred_data, Rtal.phat, Rtal.ps(1:100,:), talaris_fields);
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 12 10];
print('zemu_sigmoid2_tal_feats.png','-dpng','-r300');

plotsample(Rtal);
print('zemu_sigmoid2_tal_posterior.png','-dpng','-r600');

fit_terms = fitterms(tal_pred_data, Rtal.phat, talaris_fields, exp_data);
writetable( fit_terms, 'tal_GAM_terms.csv' )


%% control

[m,n] = size(talaris_fields)
control_pred_data = zeros( 1240, n );
for i = 1:n
    field_name = char(talaris_fields(i));
    control_pred_data(:,i) = control_table.(field_name);
end
assert( isequal( exp_data, control_table.ExperimentalDDG ) );
exp_data = control_table.ExperimentalDDG;

Rcon = sigmoidfit(control_pred_data, exp_data, 50, 1000);

plotsigmoids(control_pred_data, Rcon.phat, Rcon.ps(1:100,:), talaris_fields);
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 12 10];
print('zemu_sigmoid2_con_feats.png','-dpng','-r300');

plotsample(Rcon);
print('zemu_sigmoid2_con_posterior.png','-dpng','-r300');

fit_terms = fitterms(control_pred_data, Rcon.phat, talaris_fields, exp_data);
writetable( fit_terms, 'control_GAM_terms.csv' )



%% plot correlations from all
current_palette = { [0.29803921568627451 0.44705882352941179 0.69019607843137254], [0.33333333333333331 0.6588235294117647 0.40784313725490196], [0.7686274509803922 0.30588235294117649 0.32156862745098042], [0.50588235294117645 0.44705882352941179 0.69803921568627447], [0.80000000000000004 0.72549019607843135 0.45490196078431372] };

subplot(321);
plot(exp_data, control_table.total, 'r.', 'color', cell2mat( current_palette(2)));
title('No backrub control');
xlim([-6 11]); ylim([-6 11]);
hold on;
    coef_fit = polyfit( exp_data, control_table.total, 1 );
    y_fit = polyval( coef_fit, xlim );
    plot( xlim, y_fit, 'r', 'color', cell2mat( current_palette(2)) );
hold off;
grid on;
xlabel( 'Experimental $\Delta\Delta G$', 'Interpreter', 'latex' );
ylabel( 'Rosetta Score', 'Interpreter', 'latex' );

subplot(323);
plot(exp_data, talaris_table.total, 'r.', 'color', cell2mat( current_palette(1)));
title('Flex ddG');
xlim([-6 11]); ylim([-6 11]);
hold on;
    coef_fit = polyfit( exp_data, talaris_table.total, 1 );
    y_fit = polyval( coef_fit, xlim );
    plot( xlim, y_fit, 'r', 'color', cell2mat( current_palette(1)) );
hold off;
grid on;
xlabel( 'Experimental $\Delta\Delta G$', 'Interpreter', 'latex' );
ylabel( 'Rosetta Score', 'Interpreter', 'latex' );

subplot(325);
plot(exp_data, ref_table.total, 'r.', 'color', cell2mat( current_palette(5)));
title('Flex ddG (REF energy)');
xlim([-6 11]); ylim([-6 11]);
hold on;
    coef_fit = polyfit( exp_data, ref_table.total, 1 );
    y_fit = polyval( coef_fit, xlim );
    plot( xlim, y_fit, 'r', 'color', cell2mat( current_palette(5)) );
hold off;
grid on;
xlabel( 'Experimental $\Delta\Delta G$', 'Interpreter', 'latex' );
ylabel( 'Rosetta Score', 'Interpreter', 'latex' );

subplot(322);
plot(exp_data, Rcon.fhat, 'r.', 'color', cell2mat( current_palette(2)));
hold on;
plot( [exp_data exp_data]', [min(Rcon.fs')' max(Rcon.fs')']', '-', 'color', 0.65*[1 1 1]);
plot( exp_data, Rcon.fhat, 'r.', 'color', cell2mat( current_palette(2)));
hold off;
title('GAM No backrub control');
xlim([-6 11]); ylim([-6 11]);
hold on;
    coef_fit = polyfit( exp_data, Rcon.fhat, 1 );
    y_fit = polyval( coef_fit, xlim );
    plot( xlim, y_fit, 'r', 'color', cell2mat( current_palette(2)) );
hold off;
grid on;
xlabel( 'Experimental $\Delta\Delta G$', 'Interpreter', 'latex' );
ylabel( 'GAM Score', 'Interpreter', 'latex' );

subplot(324);
plot( exp_data, Rtal.fhat, 'r.', 'color', cell2mat( current_palette(1)));
hold on;
plot( [exp_data exp_data]', [min(Rtal.fs')' max(Rtal.fs')']', '-', 'color', 0.65*[1 1 1]);
plot( exp_data, Rtal.fhat, 'r.', 'color', cell2mat( current_palette(1)));
hold off;
title('GAM flex ddG');
xlim([-6 11]); ylim([-6 11]);
hold on;
    coef_fit = polyfit( exp_data, Rtal.fhat, 1 );
    y_fit = polyval( coef_fit, xlim );
    plot( xlim, y_fit, 'r', 'color', cell2mat( current_palette(1)) );
hold off;
grid on;
xlabel( 'Experimental $\Delta\Delta G$', 'Interpreter', 'latex' );
ylabel( 'GAM Score', 'Interpreter', 'latex' );

subplot(326);
plot( exp_data, Rref.fhat, 'r.', 'color', cell2mat( current_palette(5)));
title('GAM flex ddG (REF)');
hold on;
plot( [exp_data exp_data]', [min(Rref.fs')' max(Rref.fs')']', '-', 'color', 0.65*[1 1 1]);
plot( exp_data, Rref.fhat, 'r.', 'color', cell2mat( current_palette(5)));
hold off;
xlim([-6 11]); ylim([-6 11]);
hold on;
    coef_fit = polyfit( exp_data, Rref.fhat, 1 );
    y_fit = polyval( coef_fit, xlim );
    plot( xlim, y_fit, 'r', 'color', cell2mat( current_palette(5)) );
hold off;
grid on;
xlabel( 'Experimental $\Delta\Delta G$', 'Interpreter', 'latex' );
ylabel( 'GAM Score', 'Interpreter', 'latex' );


fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8.5 9.5];
print('zemu_sigmoid2_corrs.png','-dpng','-r600');


%% Output table of fit results


return

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
