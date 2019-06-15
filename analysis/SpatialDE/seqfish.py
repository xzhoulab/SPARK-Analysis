import pandas as pd
import NaiveDE
import SpatialDE


counts = pd.read_csv('./processed_data/seqFISH_field43_countdata.csv', index_col=0)
counts = counts.T[counts.sum(0) >= 3].T 
sample_info = pd.read_csv('./processed_data/seqFISH_field43_info.csv', index_col=0)
sample_info['total_counts'] = counts.sum(1)
counts = counts.loc[sample_info.index]  
norm_expr = NaiveDE.stabilize(counts.T).T
resid_expr = NaiveDE.regress_out(sample_info, norm_expr.T, 'np.log(total_counts)').T
sample_resid_expr=resid_expr
X = sample_info[['x', 'y']]
results = SpatialDE.run(X, sample_resid_expr)
results.to_csv('./output/seqFISH_field43_spe.csv',sep=' ', index=False, header=True)

de_results = results[(results.qval < 0.05)].copy()
ms_results = SpatialDE.model_search(X, resid_expr, de_results)
ms_results.to_csv('./output/seqFISH_field43_ms_spe.csv',sep=' ', index=False, header=True)



