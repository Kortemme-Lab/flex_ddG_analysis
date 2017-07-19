SELECT PredictionPPI.ID AS 'PredictionID', PredictionPPI.PredictionSet AS 'PredictionRunName',
UserPPDataSetExperiment.ID AS 'DataSetID', UserPPAnalysisSet.Subset,
UserPPDataSetExperiment.PDBFileID,
PredictionPPI.DDGTime AS 'AvgRuntime', PredictionPPI.maxvmem AS 'MaxMemGB',

PredictionPPIStructureScore.ScoreType, PredictionPPIStructureScore.StructureID,
PredictionPPIStructureScore.ScoreMethodID,
group_concat(PPMutagenesisMutation.RecordKey SEPARATOR ';') AS 'Mutations',
# PPMutagenesis.SKEMPI_KEY AS 'Mutations',
AVG(IFNULL(PositiveDDG.DDG, 0) - IFNULL(NegativeDDG.DDG, 0)) AS ExperimentalDDG,
PredictionPPIStructureScore.total,
(fa_atr+fa_dun+ fa_elec+ fa_intra_rep+ fa_rep+ fa_sol+ hbond_bb_sc+ hbond_lr_bb+ hbond_sc+ hbond_sr_bb+ omega+ p_aa_pp+ pro_close+ rama+ ref+ yhh_planarity) AS talaris_total_check,
(fa_atr+fa_dun+ fa_elec+ fa_intra_rep+ fa_rep+ fa_sol+ hbond_bb_sc+ hbond_lr_bb+ hbond_sc+ hbond_sr_bb+ omega+ p_aa_pp+ pro_close+ ref+ yhh_planarity+ fa_intra_sol_xover4+ rama_prepro+ lk_ball_wtd) AS ref_total_check,
fa_atr, fa_dun, fa_elec, fa_intra_rep, fa_rep, fa_sol, hbond_bb_sc, hbond_lr_bb, hbond_sc, hbond_sr_bb, omega, p_aa_pp, pro_close, rama, ref, yhh_planarity, fa_dun_dev, fa_dun_rot, fa_dun_semi, fa_intra_atr_xover4, fa_intra_elec, fa_intra_rep_xover4, fa_intra_sol_xover4, hxl_tors, lk_ball, lk_ball_bridge, lk_ball_bridge_uncpl, lk_ball_iso, rama_prepro, cart_bonded, lk_ball_wtd

FROM PredictionPPI
INNER JOIN UserPPDataSetExperiment ON PredictionPPI.UserPPDataSetExperimentID=UserPPDataSetExperiment.ID
INNER JOIN PredictionPPIStructureScore ON PredictionPPIStructureScore.PredictionPPIID=PredictionPPI.ID
INNER JOIN UserPPAnalysisSet ON UserPPAnalysisSet.UserPPDataSetExperimentID=PredictionPPI.UserPPDataSetExperimentID
INNER JOIN PPMutagenesis ON PPMutagenesis.ID=UserPPAnalysisSet.PPMutagenesisID
INNER JOIN PPMutagenesisMutation ON PPMutagenesisMutation.PPMutagenesisID=PPMutagenesis.ID
LEFT JOIN PPIDDG AS PositiveDDG ON UserPPAnalysisSet.PositiveDependentPPIDDGID=PositiveDDG.ID
LEFT JOIN PPIDDG AS NegativeDDG ON UserPPAnalysisSet.NegativeDependentPPIDDGID=NegativeDDG.ID
# INNER JOIN UserPPDataSetExperimentTag ON UserPPDataSetExperimentTag.UserPPDataSetExperimentID=UserPPDataSetExperiment.ID
# WHERE UserPPDataSetExperimentTag.Tag='ZEMu lite'
WHERE PredictionPPI.PredictionSet='%s'
# WHERE PredictionPPI.PredictionSet='zemu_1.2-60000_rscript_validated-t14'
# WHERE PredictionPPI.PredictionSet='zemu_1.2-60000_rscript_validated-ref'
# AND PredictionPPI.ID=142386 # t14 example
# AND PredictionPPI.ID=146106 # ref example
AND UserPPAnalysisSet.Subset='ZEMu'
GROUP BY PredictionPPI.ID, PredictionPPIStructureScore.StructureID, PredictionPPIStructureScore.ScoreType, PredictionPPIStructureScore.ScoreMethodID
ORDER BY PredictionPPI.ID, PredictionPPIStructureScore.StructureID, PredictionPPIStructureScore.ScoreType, PredictionPPIStructureScore.ScoreMethodID
;
