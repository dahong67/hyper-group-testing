import numpy as np
from scipy.stats import poisson
from scipy.sparse import load_npz
import glob,os
import argparse
import time
import itertools

def random_binary_balanced(m,n,q):
	A = np.zeros((m,n))
	for i in range(n):
		idx = np.random.choice(m,q,replace=False)
		A[idx,i] = 1
	return A

def random_q_pooling(m,n,q):
	s = n//(m//q)
	Atemplate = np.kron(np.eye(m//q),np.ones((1,s)))
	Alist = [Atemplate[:,np.random.permutation(n)] for i in range(q)]
	return np.vstack(Alist)

def decode(y,A,e):
	x = A.sum(0) - A.T.dot(y)
	return (x < e)

def sample_pooled_viral_load(A_i,viral_load):
	# poisson sampler has issues with super large numbers...
	viral_load[viral_load > 1e15] = 1e15
	sampled_load = poisson.rvs(viral_load/A_i.sum())
	total_sampled_load = A_i.dot(sampled_load)
	return total_sampled_load

def test_n_samples(A,viral_load,LOD,fp_rate):
	y = np.zeros(A.shape[0])
	for i in range(len(y)):
		total_sampled_load = sample_pooled_viral_load(A[i],viral_load)
		if (total_sampled_load > LOD) or (np.random.random() < fp_rate): # call the pool positive with probability fp_rate, regardless of viral load
			y[i] = 1
	return y

def run_test(A,ViralLoad,InfectionTime,OnsetTime,LOD,err_tol,fp_rate,numpos,itr):
	nperms = sum(1 for _ in itertools.permutations(range(A.shape[1]),numpos))
	sample_true_pos = np.zeros((A.shape[1],nperms))
	sample_total_pos = np.zeros((A.shape[1],nperms))
	num_stage_two = np.zeros(nperms)
	timingstart = time.time()
	for __ in range(itr):
		if __ % max(1,itr//100) == 0:
			print('Iteration: %d; Runtime: %.5f; total_pos: %d' % (__,time.time()-timingstart,sample_total_pos.sum()))
		pos_viral_load = np.random.choice(ViralLoad[ViralLoad != 0],numpos)
		for permidx,perm in enumerate(itertools.permutations(range(A.shape[1]),numpos)):
			viral_load = np.zeros(A.shape[1])
			viral_load[np.asarray(perm)] = pos_viral_load
			x = (viral_load > 0).astype(np.float)
			y = test_n_samples(A,viral_load,LOD,fp_rate)
			x_hat = decode(y,A,err_tol)
			t = x_hat.sum() # all putative positives
			num_stage_two[permidx] += t
			x_hat = x_hat*(viral_load > LOD) # samples that are putative positive and above the LOD
			sample_true_pos[:,permidx] += x*x_hat
			sample_total_pos[:,permidx] += x
	num_stage_two /= itr
	return sample_true_pos,sample_total_pos,num_stage_two

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--viral-load-matrix', help='Path to viral load matrix (individuals x time points) in sparse format')
	parser.add_argument('--infection-time', help='Path to time of infection for each individual')
	parser.add_argument('--savepath', help='Path to save summary figure')
	parser.add_argument('--n-individuals', help='Number of individuals to test (n)',type=int)
	parser.add_argument('--m-pools', help='Number of pools (m)',type=int)
	parser.add_argument('--q-split', help='Number of pools per sample (q)',type=int)
	parser.add_argument('--numpos', help='Number of positives per batch',type=int)
	parser.add_argument('--numiters', help='Number of iterations',type=int,default=1000)
	parser.add_argument('--randseed', help='Random number seed',type=int,default=0)
	parser.add_argument('--at-time', help='Time point to analyze',type=int)
	parser.add_argument('--fp-rate', help='Rate of false positive PCRs',type=float,default=0.01)
	parser.add_argument('--LOD', help='Limit of detection',type=float,default=100)
	parser.add_argument('--pool-compositions', help='Path to matrix of pool compositions, otherwise random (default)',default=None)
	parser.add_argument('--q-pooling', help='Flag to use random q-pooling, otherwise random binary balanced (default)',default=None)
	parser.add_argument('--designseed', help='Seed for random designs',type=int)
	parser.add_argument('--expected-pos-prob', help='1-expected faulty test rate (for setting error correction threshold)',type=float,default=0.95)
	args,_ = parser.parse_known_args()
	# number of faulty tests to tolerate
	err_tol = np.round(2*(1-args.expected_pos_prob)*args.q_split) + 1e-7
	for key,value in vars(args).items():
		print('%s\t%s' % (key,str(value)))
	print('error_tol\t%d' % err_tol)
	f = open(args.infection_time)
	header = f.readline()
	InfectionTime = []
	for line in f:
		if '.' in line: # non-empty lines
			InfectionTime.append(float(line.strip().split(',')[-1]))
		else:
			InfectionTime.append(-1)
	f.close()
	f.close()
	InfectionTime = np.array(InfectionTime)
	print('parsed infection times for %d individuals' % len(InfectionTime))
	ViralLoad = load_npz(args.viral_load_matrix)
	timepoints = np.load(args.viral_load_matrix.replace('viral_loads.npz', 'timepoints.npy'))
	PeakTime = np.load(args.viral_load_matrix.replace('viral_loads.npz', 'peak_times.npy'))
	print('loaded viral load matrix with shape (%d, %d) and %d entries' % (ViralLoad.shape[0], ViralLoad.shape[1], len(ViralLoad.data)))
	# pooling composition matrix
	if args.pool_compositions is not None:
		print('loading pools from: %s' % args.pool_compositions)
		A = np.load(args.pool_compositions)
	else:
		if args.q_pooling is not None:
			if (args.m_pools % args.q_split) == 0 and (args.n_individuals % (args.m_pools / args.q_split)) == 0:
				print('using random q-pooling')
				np.random.seed(args.designseed)
				A = random_q_pooling(args.m_pools,args.n_individuals,args.q_split)
			else:
				print('invalid params for random q-pooling. exiting.')
				exit()
		else:
			print('using random binary balanced')
			np.random.seed(args.designseed)
			A = random_binary_balanced(args.m_pools,args.n_individuals,args.q_split)
	t = args.at_time
	viralload = ViralLoad[:,t].toarray()[:,0]
	np.random.seed(args.randseed)
	sample_true_pos,sample_total_pos,num_stage_two = run_test(A,viralload,t-InfectionTime,t-PeakTime,args.LOD,err_tol,args.fp_rate,args.numpos,args.numiters)
	np.save('%s/sens_avg_sample.n-%d_m-%d_q-%d.t-%d.numpos-%d.numiters-%d.randseed-%d.npy' % (args.savepath,args.n_individuals,args.m_pools,args.q_split,t,args.numpos,args.numiters,args.randseed),sample_true_pos.sum(1)/sample_total_pos.sum(1))
	np.save('%s/sens_avg_perm.n-%d_m-%d_q-%d.t-%d.numpos-%d.numiters-%d.randseed-%d.npy' % (args.savepath,args.n_individuals,args.m_pools,args.q_split,t,args.numpos,args.numiters,args.randseed),sample_true_pos.sum(0)/sample_total_pos.sum(0))
	np.save('%s/num_stage_two.n-%d_m-%d_q-%d.t-%d.numpos-%d.numiters-%d.randseed-%d.npy' % (args.savepath,args.n_individuals,args.m_pools,args.q_split,t,args.numpos,args.numiters,args.randseed),num_stage_two)

