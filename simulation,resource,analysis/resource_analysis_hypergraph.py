import numpy as np
from scipy.sparse import load_npz
import glob,os
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
mpl.style.use('ggplot')
import seaborn as sns
import argparse

def bin_results(R,P0,n=50,y0=None):
	P = np.copy(P0)
	lower,upper = np.percentile(P,1),np.percentile(P,99)
	P[P<lower] = lower
	P[P>upper] = upper
	xs = np.argsort(P)
	bins = [P[xs[i]] for i in range(0,len(P),int(len(P)/n))]
	indices = np.digitize(P,bins)
	x = np.array([np.average(P[indices == i]) for i in set(indices)])
	if y0 is not None:
		y = np.array([np.average(R[indices == i] < y0) for i in set(indices)])
	else:
		y = np.array([np.average(R[indices == i]) for i in set(indices)])
	return x,y

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--viral-load-matrix', help='Path to viral load matrix (individuals x time points)')
	parser.add_argument('--resultspath', help='Path to simulation results')
	parser.add_argument('--savepath', help='Path to directory to save summary file')
	parser.add_argument('--start-time', help='Beginning of time range to analyze',type=int)
	parser.add_argument('--end-time', help='End of time range to analyze',type=int)
	parser.add_argument('--LOD', help='Limit of detection',type=float,default=100)
	parser.add_argument('--n-swabs-list',nargs='+',help='List of swab budgets',type=int)
	parser.add_argument('--m-kits-list',nargs='+',help='List of test kit budgets',type=int)
	args,_ = parser.parse_known_args()
	ViralLoad = load_npz(args.viral_load_matrix)
	timepoints = np.load(args.viral_load_matrix.replace('viral_loads.npz', 'timepoints.npy'))
	e0 = np.average([np.average(v.data > args.LOD) for v in ViralLoad[:,args.start_time:(args.end_time+1)].T if len(v.data)>0])
	Results = []
	FP = glob.glob(os.path.join(args.resultspath, '*', 'Recall_combined.*.npy'))
	for fp in FP:
		cond = fp.split('.')[-2]
		n = int(cond.split('_')[0].split('-')[1])
		m = int(cond.split('_')[1].split('-')[1])
		q = int(cond.split('_')[2].split('-')[1])
		q = min(m,q) # filenames sometime screwy
		method = fp.split('/')[-2].replace(',','_')
		sensitivity = np.load(fp)
		efficiency = np.load(fp.replace('Recall_combined','Eff_avg'))
		avg_sensitivity = sensitivity[args.start_time:(args.end_time+1)]
		avg_efficiency = efficiency[args.start_time:(args.end_time+1)]
		Results.append((n,m,q,avg_sensitivity,avg_efficiency,method))
	f = open('%s/summary.resource_t0-%d_t1-%d.csv' % (args.savepath, args.start_time, args.end_time),'w')
	_ = f.write(','.join(['Sample budget','Test budget','Design samples','Design pools','Design sample split','Design runs per day','Design method','Design effectiveness']) + '\n')
	for i,n_swabs in enumerate(args.n_swabs_list):
		for j,m_kits in enumerate(args.m_kits_list):
			nm_min = min(n_swabs,m_kits)
			best = (e0*nm_min,nm_min,nm_min,1,1,'individual')
			_ = f.write('%d,%d,%d,%d,%d,%f,%s,%f\n' % (n_swabs,m_kits,best[1],best[2],best[3],best[4],best[5],best[0]))
			best = (0,0,0,0,0,'nothing')
			for n,m,q,sens,eff,method in Results:
				m_eff = n/eff # stage (i) pools + stage (ii) validation tests
				n_tests = np.minimum(n_swabs/n, m_kits/m_eff)
				if np.average(n_tests) > 0.9:
					effective_tests = np.average(n*n_tests*sens)
					if method != 'hypergraph':
						_ = f.write('%d,%d,%d,%d,%d,%f,%s,%f\n' % (n_swabs,m_kits,n,m,q,np.average(n_tests),method,effective_tests))
					elif effective_tests > best[0]:
						best = (effective_tests,n,m,q,np.average(n_tests),method)
			if best[5] != 'nothing':
				_ = f.write('%d,%d,%d,%d,%d,%f,%s,%f\n' % (n_swabs,m_kits,best[1],best[2],best[3],best[4],best[5],best[0]))
	f.close()

# To run:
# conda activate covid-group-tests
# python -u resource_analysis_hypergraph.py --viral-load-matrix ../simulation/Simulated_populations/seir_viral_loads_swab.viral_loads.npz --resultspath ../simulation/designs/ --savepath . --start-time 40 --end-time 90 --n-swabs-list 12 24 48 96 192 384 768 1536 3072 6144 --m-kits-list 12 24 48 96 192 384 768 1536 3072 6144