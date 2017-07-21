SELECT UserPPDataSetExperiment.ID AS 'DataSetID', UserPPAnalysisSet.Subset,
UserPPDataSetExperiment.PDBFileID,
PDBFile.Resolution,

group_concat(PPMutagenesisMutation.RecordKey SEPARATOR ';') AS 'Mutations',
AVG(IFNULL(PositiveDDG.DDG, 0) - IFNULL(NegativeDDG.DDG, 0)) AS ExperimentalDDG

FROM UserPPDataSetExperiment
INNER JOIN UserPPAnalysisSet ON UserPPAnalysisSet.UserPPDataSetExperimentID=UserPPDataSetExperiment.ID
INNER JOIN PPMutagenesis ON PPMutagenesis.ID=UserPPAnalysisSet.PPMutagenesisID
INNER JOIN PPMutagenesisMutation ON PPMutagenesisMutation.PPMutagenesisID=PPMutagenesis.ID
INNER JOIN PDBFile ON PDBFile.ID=UserPPDataSetExperiment.PDBFileID
LEFT JOIN PPIDDG AS PositiveDDG ON UserPPAnalysisSet.PositiveDependentPPIDDGID=PositiveDDG.ID
LEFT JOIN PPIDDG AS NegativeDDG ON UserPPAnalysisSet.NegativeDependentPPIDDGID=NegativeDDG.ID

WHERE UserPPAnalysisSet.Subset='ZEMu'
GROUP BY DataSetID
ORDER BY DataSetID
;
