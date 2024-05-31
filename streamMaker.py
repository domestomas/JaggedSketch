#!/usr/bin/python3
# Licensed to the Apache Software Foundation (ASF) under one or more
# contributor license agreements.  See the NOTICE file distributed with
# this work for additional information regarding copyright ownership.
# The ASF licenses this file to You under the Apache License, Version 2.0
# (the "License"); you may not use this file except in compliance with
# the License.  You may obtain a copy of the License at
#
#	 http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

'''
Written by Edo Liberty and Pavel Vesely.
'''

import random
from math import sqrt,ceil

class StreamMaker():
	def __init__(self):
		self.orders = ['sorted','reversed','zoomin','zoomout','sqrt','random','adv','clustered', 'clustered-zoomin', 'random'] 
		
	def make(self, n=1000, order='random', p=1000, g=0, s=1):
		#assert(order in self.orders)
		assert order in self.orders
		
		if order == 'sorted': # sorted order
			for item in range(1, n+1):
				yield item
		elif order == 'reversed': # reversed sorted order
			for item in range(n, 0, -1):
				yield item
		elif order == 'zoomin': 
			for item in range(1, int(n/2+1)):
				yield item
				yield n-item+1
		elif order == 'zoomout': 
			for item in range(0, int(n/2)):
				yield n/2 + item+1
				yield n/2 - item
		elif order == 'sqrt': 
			t = int(sqrt(2*n))
			item = 0
			initialItem = 0
			initialSkip = 1
			for i in range(t):
				item = initialItem
				skip = initialSkip
				for j in range(t-i):
					yield item+1
					item+=skip
					skip+=1 
				initialSkip+=1
				initialItem+=initialSkip
		elif order == 'adv':
			m = ceil(n/p)
			for i in range(p):
				for j in range(s*(g + p + m*(p-i)), s*(g + p + m*(p-i-1)), -s):
					yield j
				yield i
				if i == p // 2:
					for j in range(p, s*(g + p + m), s*(g + p + m) // 10):
						yield j
		elif order == 'clustered': # sorted clustered order
			m = ceil(n/p) # number of clusters  of size p
			for i in range(m):
				# output cluster; g is the gap between clusters
				for j in range(i*g, i*g + p):
					yield i*g + j / p
			for i in range(m):
				# put some items (roughly s many) into the gap between clusters
				for j in range(i*g + p, (i+1)*g, g // s):
					yield j
		elif order == 'clustered-zoomin': # order roughly as in zoomin
			m = ceil(n/p) # number of clusters  of size p
			for i in range(m):
				# output cluster; g is the gap between clusters
				for j in range(i*g, i*g + p, 2):
					yield i*g + j / p
			for i in range(m):
				# put some items (roughly s many) into the gap between clusters
				for j in range(i*g + p, (i+1)*g, g // s):
					yield j
			for i in range(m - 1, 0, -1):
				# output cluster; g is the gap between clusters
				for j in range(i*g + p, i*g, -2):
					yield i*g + (j + 1) / p
		else: # order == 'random':
			for item in random.sample(range(1, n+1), k=n):
				yield item

if __name__ == '__main__':
	import sys
	import argparse
	streamer = StreamMaker()

	parser = argparse.ArgumentParser()
	parser.add_argument(
		'-n', type=int, default=1000,
		help='the number of generated elements'
	)
	parser.add_argument(
		'-p', type=int, default=1000,
		help='parameter for generating some orders (for orders adv and clustered only)'
	)
	parser.add_argument(
		'-g', type=int,  default=0,
		help='another parameter for generating some orders (for orders adv and clustered only)'
	)
	parser.add_argument(
		'-s', type=int, default=1,
		help='yet another parameter for generating some orders (for orders adv and clustered only)'
	)
	parser.add_argument(
		'-order', type=str, choices=streamer.orders, default="random",
		help='the order of the streamed integers.',
	)
	args = parser.parse_args()
	
	for item in streamer.make(args.n, args.order, args.p, args.g, args.s):
		sys.stdout.write('%d\n'%item)