#!/usr/bin/env python3
from streamMaker import StreamMaker
from jaggedSketchSimple import JaggedSketch
import argparse
from collections import namedtuple
import multiprocessing as mp
from multiprocessing import shared_memory
import pickle, os, random
import numpy as np
from pprint import pprint
from random import randrange


Data = namedtuple('Data', ['perc68', 'perc95', 'perc99', 'avg', 'median'])
def run_the_sketch(n, order, q, J, epsilon, mem_name=''):
		sketch = JaggedSketch(epsilon=epsilon, important_quantiles=q, constant_J = J)
		if mem_name != '':
			mem = shared_memory.SharedMemory(mem_name)
			stream = np.ndarray((n,), dtype=np.int64, buffer=mem.buf)
		else:
			stream = StreamMaker().make(n=n, order=order)
		
		for item in stream:
			sketch.update(item)
			
		if mem_name != '':
			mem.close()
		return sketch

def bisect(n, order, q, j, space):
	cap = 0
	small = 0.001
	big = 0.1
	while big-small > 0.00001 and abs(cap-space) > 10:
		avg = round((small + big) / 2, 6)
		s = run_the_sketch(n, 'sorted', q, j, avg)
		cap = sum(c.capacity for c in s.compactors)
		if cap > space:
			small = avg
		else:
			big = avg
		#print(f"cap: {cap} | avg: {avg} small: {small} big: {big}")
	return avg

def main():
	parser = argparse.ArgumentParser(description=
		'Program for testing Jagged Sketch.'
	)
	parser.add_argument(
		'-n', type=int, default=10000000, 
		help='the approximate number of generated elements'
	)
	parser.add_argument(
		'-order', type=str, default='random',
		choices = ['sorted','reversed','zoomin','zoomout','sqrt', 'random'], 
		help='the order of the streamed integers (in -random case stream must fit in memory)'
	)
	parser.add_argument(
		'-repeat', type=int, default=1000,
		help='the number of times to repeat building the sketch (to calculate deviations)'
	)
	parser.add_argument(
		'-info', type=str, default='',
		help='additional info to include in the legend'
	)
	parser.add_argument(
		'-q', type=float, default=[], action='append',
		help='important quantiles, default = {}'
	)
	parser.add_argument(
		'-j', type=float, default=0.5,
		help='the constant J from theory'
	)
	parser.add_argument(
		'-epsilon', type=float, default=0, 
	)
	parser.add_argument(
		'-space', type=float, default=10020, 
	)
	args = parser.parse_args()
	
	# parse input
	n = args.n
	repeat = args.repeat
	order = args.order
	user_info = f"{args.info+'_' if args.info!='' else ''}simple"
	epsilon = args.epsilon
	space = args.space
	q = set(args.q) if args.q != [] else {0}
	J = args.j
	
	if epsilon == 0:
		epsilon = bisect(n, order, q, J, space)
	
	filename = (
		f"js_{int(n/1000000)}mil_{order}"
		f"_q{''.join([str(x) for x in q])}_J{J}"
		f"_eps{epsilon}{'_'+user_info if user_info!= '' else ''}"
		)
	
	# check files and folders
	os.makedirs("samples", exist_ok=True)
	if repeat > 1 and os.path.isfile(f"samples/{filename}"):
		exit("file already exists")
	
	mem_name = ''
	if order == "random":
		mem_name = 'my_mem_'+str(randrange(10000))
		rng = np.random.default_rng()
		a = np.array(rng.permutation(np.arange(1,n+1)))
		mem = shared_memory.SharedMemory(create=True, size=a.nbytes, name=mem_name)
		random_stream = np.ndarray((n,), np.int64, buffer=mem.buf)
		random_stream[:] = a[:]
	
	# run the sketch in paralel "repeat" times
	with mp.Pool() as pool:
		async_runs = [pool.apply_async(
				run_the_sketch, (n, order, q, J, epsilon, mem_name)
			) for _ in range(repeat)]
		runs = [x.get() for x in async_runs]
	
	if order == "random":
		mem.close()
		mem.unlink()
	
	# get skech info
	s = runs[0]
	n = s.N
	sketch_info = {
		"n":s.N, "repeat":repeat, "cap":sum(c.capacity for c in s.compactors), 
		"B":max(c.capacity for c in s.compactors), "H":s.H(), "J":s.J, 
		"epsilon":s.epsilon, "Q":s.important_quantiles, "user":user_info
	}
	
	if repeat == 1:
		pprint(sketch_info)
		for c in s.compactors:
			print(f"ss: {c.section_size}, comp: {c.num_compactions} B: {c.capacity}")

	ranks = [r.ranks() for r in runs]
	
	# dump the results to file
	if repeat > 1:
		with open(f"samples/{filename}", mode='xb') as file:
			pickle.dump(Sampling(ranks, sketch_info, n), file)

class Sampling():

	def __init__(self, ranks, info, n):
		self.info = info
		self.repeat = len(ranks)
		self.sample_points = self.choose_sample_points(ranks)
		self.data = Data([], [], [], [], [])
		self.n = n
		self.prepare_data(ranks)
	
	def choose_sample_points(self, ranks):
		all_points = [x for run in ranks for (x, _) in run]
		all_points.sort()
		return [all_points[i*self.repeat+1] for i in range(len(all_points)//(self.repeat)-1)]

	def prepare_data(self, ranks):
		repeat = self.repeat
		current_index = [0 for _ in range(repeat)]
		sd_percents = (68*repeat//100, 95*repeat//100, 99*repeat//100)
		
		# find SD values for each point
		# go through all the lists in sync to achieve linear time
		for point in self.sample_points:
			errors = []
			error_sum = 0
			
			# for every run
			for i in range(repeat):
				run = ranks[i]
				curr = current_index[i]
				
				# move a little to find error in given point
				while curr < len(run) and run[curr][0] <= point: 
					curr += 1
				current_index[i] = curr
				
				# calculate error in given point
				err = run[curr-1][1] - point

				# save the error
				error_sum += err
				errors.append(err)
			
			# find an save the confidence intervals
			errors.sort()
			self.data.median.append(errors[len(errors)//2])
			self.data.avg.append(error_sum/repeat)
			errors = [abs(e) for e in errors]
			errors.sort()
			for i in range(3):
				self.data[i].append(errors[sd_percents[i]])

if __name__ == '__main__':
	main()
