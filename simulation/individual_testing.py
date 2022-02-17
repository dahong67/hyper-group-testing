import numpy as np
from scipy.sparse import load_npz
import argparse

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--viral-load-matrix', help='Path to viral load matrix (individuals x time points) in sparse format')
	parser.add_argument('--savepath', help='Path to save summary figure')
	parser.add_argument('--start-time', help='Start of time range to analyze',type=int)
	parser.add_argument('--end-time', help='End of time range to analyze',type=int)
	parser.add_argument('--LOD', help='Limit of detection',type=float,default=100)
	args,_ = parser.parse_known_args()
	# Load data
	for key,value in vars(args).items():
		print('%s\t%s' % (key,str(value)))
	ViralLoad = load_npz(args.viral_load_matrix)
	print('loaded viral load matrix with shape (%d, %d) and %d entries' % (ViralLoad.shape[0], ViralLoad.shape[1], len(ViralLoad.data)))
	# Compute prevalence, average efficiency and sensitivity
	Prevalence = np.zeros(ViralLoad.shape[1])
	E_avg = np.zeros(ViralLoad.shape[1])
	Recall_combined = np.zeros(ViralLoad.shape[1])
	for t in range(args.start_time,args.end_time+1):
		viralload = ViralLoad[:,t].toarray()[:,0]
		x_true = (viralload > 0).astype(np.float)
		x_hat = (viralload > args.LOD).astype(np.float)
		prev = np.average(viralload > 0)
		ea = 1.0
		true_pos = x_true.dot(x_hat)
		total_pos = x_true.sum()
		if total_pos > 0:
			rec = true_pos/total_pos
		else:
			rec = 1
		Prevalence[t] = prev
		E_avg[t] = ea
		Recall_combined[t] = rec
		print('Time: %d; Case frequency: %.5f; Efficiency: %.2f; Sensitivity: %.4f' % (t,prev,ea,rec))
	np.save('%s/Eff_avg.npy' % (args.savepath),E_avg)
	np.save('%s/Recall_combined.npy' % (args.savepath),Recall_combined)
	np.save('%s/Prevalence.npy' % (args.savepath),Prevalence)

# python -u individual_testing.py --viral-load-matrix Simulated_populations/seir_viral_loads_swab.viral_loads.npz --savepath designs/individual --start-time 0 --end-time 356 | tee designs/individual/log.txt
