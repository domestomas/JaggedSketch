#!/usr/bin/python3

from random import random
from math import log

# CONSTANTS
SMALLEST_MEANINGFUL_SECTION_SIZE = 4
INIT_SECTIONS = 1.5

class JaggedSketch:
	def __init__(self, epsilon=0.01, delta=0.01, important_quantiles={0}, 
			constant_J=0.5, improvement_for_high_ranks=True):
		if epsilon <= 0 or epsilon > 1:
			raise ValueError("epsilon must be between 0 and 1")
		if delta <= 0 or delta > 0.5:
			raise ValueError("delta must be between 0 and 0.5")
		if constant_J < 0:
			raise ValueError("J must be non-negative")
		if not all(x >= 0 and x <= 1 for x in important_quantiles):
			raise ValueError("All important quantiles must be between 0 and 1")
		if constant_J != 0 and important_quantiles == set():
			raise ValueError("with no important quentiles, j must equal 0")
		# Set of ranks with higher accuracy given as an imput to the quantile function
		self.important_quantiles = important_quantiles
		# Gives importance of ranks in Q; 
  		# J=0 means all ranks have the same importance
		self.J = constant_J
		# Relative error for the desired rank in Q (quarantee from the theory)
		self.epsilon = epsilon
		# Error improvement for high ranks
		self.improvement_for_high_ranks = improvement_for_high_ranks
		# Delta is the probability of error larger than epsilon for given query
		self.probability_constant = log(1/delta)**0.5		
		# Size of the input summarized
		self.N = 0
		# Current number of saved items
		self.size = 0
		# Sum of capacities of all compactors
		self.capacity = 0
		# Levels corresponding to important quantiles
		self.important_levels = set()
		self.compactors = []
		self.compactors.append(RelativeCompactor(self))
		self.compactors[0].set_capacity_and_section_size()
	
	def H(self):
		return len(self.compactors)

	# Adds new compactor to the sketch
	def grow(self):
		# Add a new compactor
		self.compactors.append(RelativeCompactor(self))
		self.compactors[-1].set_capacity_and_section_size()
		
		# Do the full compaction for all compactors
		for (h, compactor) in enumerate(self.compactors[:self.H()-1]):
			self.compactors[h+1].extend(compactor.full_compaction())
		while self.compactors[-1].is_full():
			self.compactors.append(RelativeCompactor(self))
			self.compactors[-1].set_capacity_and_section_size()
			self.compactors[-1].extend(self.compactors[-2].full_compaction())
		
		# Update all the parameters
		self.update_important_levels()
		for c in self.compactors:
			c.set_capacity_and_section_size()
			
	# Adds new item to the skech
	def update(self, item):
		self.compactors[0].append(item)
		self.N += 1
		self.size += 1
		if self.size >= self.capacity:
			self.compress()
		assert self.size < self.capacity
	
	# Do the compaction on level zero and possibly on higher levels
	def compress(self):
		for (h, compactor) in enumerate(self.compactors):
			if compactor.is_full():
				if h+1 == self.H():
					self.grow()
					return
				self.compactors[h+1].extend(compactor.normal_compaction())
				# Be lazy and do not continue under capacity
				if self.size < self.capacity:
					return
				
	
	# Find the right levels corresponding to the quantiles
	# We assume that when this function is called, all compactors are sorted
	def update_important_levels(self):
		self.important_levels.clear()
		for q in self.important_quantiles:
			x = self.quantile(q) # item with appropriate quantile
			
			# binary seach for the right level
			i = 0
			j = self.H() - 1
			while i < j-1:
				m = (i + j) // 2
				level_min = self.compactors[m][0]
				if x >= level_min:
					i = m
				else:
					j = m

			# save the calculated level
			self.important_levels.add(i)
	
	# Computes a list of items and their ranks
	def ranks(self):
		ranks_list = []
		items_and_weights = []
		for (h, items) in enumerate(self.compactors):
			items_and_weights.extend( (item, 2**h) for item in items )
		items_and_weights.sort()
		cum_weight = 0
		for (item, weight) in items_and_weights:
			cum_weight += weight
			ranks_list.append( (item, cum_weight) )
		return ranks_list

	# Computes cummulative distribution function (as a list of items 
	# and their ranks expressed as a number in [0,1])
	def cdf(self):
		cdf = []
		rank_list = self.ranks()
		_, total_weight = rank_list[-1]
		for (item, cum_weight) in rank_list:
			cdf.append( (item, cum_weight / total_weight) )
		return cdf

	# Returns an approximate rank of value
	def rank(self, value):
		return sum(c.rank(value)*2**h for (h, c) in enumerate(self.compactors))
		

	# Returns an input item which is approx. q-quantile 
 	# (i.e. has rank approx. q*self.N)
	def quantile(self, q):
		assert (q >= 0 and q <= 1), f"parameter q must be in [0, 1], but q = {q}"
		desired_rank = q*self.N
		ranks = self.ranks()
		i = 0
		j = len(ranks)
		while i < j:
			m = (i + j) // 2
			(item, rank) = ranks[m]
			if desired_rank > rank:
				i = m + 1
			else: j = m
		(item, rank) = ranks[i]
		return item

class RelativeCompactor(list):
	def __init__(self, sketch):
		self.num_compactions = 0 # Number of compaction operations performed
		self.state = 0 # State of the deterministic compaction schedule
		self.offset = 0 # Indicator for taking even or odd items
		self.shift = 0 # Indicator for shifting the compacted part by one item
		self.sketch = sketch
		self.h = sketch.H() # height (level) of the compactor
		self.capacity = None
		self.section_size = None
		
	def rank(self, value):
		return sum(1 for v in self if v <= value)

	def is_full(self):
		return len(self) >= self.capacity
	
	def reset_compaction_schedule(self):
		self.state = 0
		self.set_capacity_and_section_size()
	
	def set_capacity(self):
		old_capacity = self.capacity if self.capacity != None else 0
		if self.sketch.improvement_for_high_ranks:
			self.capacity = int(self.sketch.probability_constant * 
				self.sketch.H()**(0.5 + min(1, self.sketch.J)) /
				(self.scale() * self.sketch.epsilon)
			)
		else:
			self.capacity = int(self.sketch.probability_constant * 
				self.sketch.H()**min(1, self.sketch.J) * 
				log(2 + self.num_compactions, 2)**0.5 /
				(self.scale() * self.sketch.epsilon)
			)
		self.sketch.capacity += self.capacity - old_capacity
		
	def set_section_size(self):
		self.section_size = int(
			self.capacity / 
			( 2 * INIT_SECTIONS * log(2 + self.num_compactions, 2) )
		)
	
	def set_capacity_and_section_size(self):
		self.set_capacity()
		self.set_section_size()

	# Chooses a scaling factor by the distance to the closest important level
	def scale(self):
		# distance from the closest important level
		if len(self.sketch.important_levels) > 0:
			dist = min([abs(self.h-l) for l in self.sketch.important_levels])
		else:
			dist = 0
		
		# choose the scaling factor (based on J parameter)
		if dist == 0: 
			scale = 1
		elif dist == 1: 
			scale = 1.5**self.sketch.J
		else: 
			scale = dist**self.sketch.J
		
		return min(scale, self.sketch.H())
	
	# Counts the number of protected items based on the compaction schedule
	def count_protected(self):
		right_part = self.capacity // 2
		rest = len(self) - self.capacity
		section_size = self.section_size
		
		# If the section size is too small we do not use the schedule
		if section_size < SMALLEST_MEANINGFUL_SECTION_SIZE: 
			compacted = right_part + rest
		else: 
			sections_to_compact = trailing_ones_binary(self.state) + 1
			self.state += 1
			right_compacted = sections_to_compact * section_size
			# schedule overflow
			if right_compacted >= right_part:
				right_compacted = right_part
				self.reset_compaction_schedule()
			compacted = right_compacted + rest
		compacted += compacted % 2
		return len(self) - compacted
	
	# Compacts everything except for the left half and resets the schedule
	def full_compaction(self):
		protected = self.capacity // 2 + 1
		protected -= (len(self)-protected) % 2
		self.reset_compaction_schedule()
		for x in self.compact(protected):
			yield x

	# Standard compaction by the schedule	
	def normal_compaction(self):
		assert self.is_full()
		for x in self.compact(self.count_protected()):
			yield x
	
	# Compacts all items exept the smallest "protected"
	def compact(self, protected):
		assert len(self[protected: ]) % 2 == 0
		self.sort() 
		
		# Set the random offset and random shift independently
		# each choice every other time
		if self.num_compactions % 2 == 1:
			self.offset = 1 - self.offset
			self.shift = int(random() < 0.5)
		else:
			self.offset = int(random() < 0.5)
			self.shift = 1 - self.shift
		
		# yield half of non-protected and delete all of them from self
		for i in range(protected + self.offset - self.shift, len(self) - self.shift, 2):
			yield self[i] # yield selected items
		self.sketch.size -= len(self[protected: ]) // 2
		del self[protected - self.shift : len(self) - self.shift]
		self.num_compactions += 1
		assert not self.is_full()

# AUXILIARY FUNCTIONS
def trailing_ones_binary(n):
	s = str("{0:b}".format(n))
	return len(s)-len(s.rstrip('1'))