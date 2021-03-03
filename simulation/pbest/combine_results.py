import numpy as np
import argparse

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--resultspath', help='Path to simulation results')
	parser.add_argument('--savepath', help='Path to directory to save combined results file')
	parser.add_argument('--start-time', help='Beginning of time range to analyze',type=int)
	parser.add_argument('--end-time', help='End of time range to analyze',type=int)
	args,_ = parser.parse_known_args()
	n = 384
	m = 48
	q = 6
	
	E_avg = np.zeros(args.end_time+1)
	Recall_combined = np.zeros(args.end_time+1)
	Specificity_combined = np.zeros(args.end_time+1)
	for t in range(args.start_time,args.end_time+1):
		ea = np.load('%s/Eff_avg.n-%d_m-%d_q-%d.t-%d.npy' % (args.resultspath,n,m,q,t)).item()
		rec = np.load('%s/Recall_combined.n-%d_m-%d_q-%d.t-%d.npy' % (args.resultspath,n,m,q,t)).item()
		spec = np.load('%s/Specificity_combined.n-%d_m-%d_q-%d.t-%d.npy' % (args.resultspath,n,m,q,t)).item()
		E_avg[t] = ea
		Recall_combined[t] = rec
		Specificity_combined[t] = spec

	np.save('%s/Eff_avg.n-%d_m-%d_q-%d.npy' % (args.savepath,n,m,q),np.array(E_avg))
	np.save('%s/Recall_combined.n-%d_m-%d_q-%d.npy' % (args.savepath,n,m,q),np.array(Recall_combined))
	np.save('%s/Specificity_combined.n-%d_m-%d_q-%d.npy' % (args.savepath,n,m,q),np.array(Specificity_combined))
