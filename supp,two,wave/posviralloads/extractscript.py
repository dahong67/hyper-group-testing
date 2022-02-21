import numpy as np
from scipy.sparse import load_npz
import argparse

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--viral-load-matrix', help='Path to viral load matrix (individuals x time points) in sparse format')
	parser.add_argument('--savepath', help='Path to save positive viral loads')
	parser.add_argument('--at-time', help='Time point to analyze',type=int)
	args,_ = parser.parse_known_args()
	for key,value in vars(args).items():
		print('%s\t%s' % (key,str(value)))
	ViralLoad = load_npz(args.viral_load_matrix)
	print('loaded viral load matrix with shape (%d, %d) and %d entries' % (ViralLoad.shape[0], ViralLoad.shape[1], len(ViralLoad.data)))
	# 
	t = args.at_time
	viralload = ViralLoad[:,t].toarray()[:,0]
	posviralload = viralload[viralload.nonzero()]
	np.save('%s/posviralloads.t-%d.npy' % (args.savepath,t), posviralload)

# python -u extractscript.py --viral-load-matrix ../Simulated_populations_two_wave/seir_viral_loads_swab_switch_SEIR_new_1.viral_loads.npz --savepath . --at-time 82
# python -u extractscript.py --viral-load-matrix ../Simulated_populations_two_wave/seir_viral_loads_swab_switch_SEIR_new_1.viral_loads.npz --savepath . --at-time 104
# python -u extractscript.py --viral-load-matrix ../Simulated_populations_two_wave/seir_viral_loads_swab_switch_SEIR_new_1.viral_loads.npz --savepath . --at-time 173