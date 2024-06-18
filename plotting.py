#!/usr/bin/env python3
import matplotlib
import matplotlib.pyplot as plt
from pprint import pformat
import sys, os, pickle, subprocess
from JSSTest import Sampling, Data

def main():
	if sys.argv[1] == "-":
		files = sys.stdin.read().rstrip().splitlines()
	else:
		files = sys.argv[1:]
	Plotting.plot(files)

class Plotting():
	def plot(files):
		rel_fig, rel_ax = plt.subplots()
		add_fig, add_ax = plt.subplots()
		colours = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'grey', 'orange', 'purple', 'brown', 'lightgreen']
		
		# plot all the datasets
		for (i, file) in enumerate(files):
			Plotting.add_dataset(file, rel_ax, add_ax, colours[i])
		
		add_ax.set_xscale('log')
		rel_ax.set_xscale('log')
		
		#
		# RELATIVE VERSION
		#
		
		# set labels and lims
		rel_ax.set_xticks([10**x for x in range(10)])
		rel_ax.set_xlabel("Rank", fontsize=10)
		rel_ax.xaxis.set_label_coords(.5, .04)
		rel_ax.yaxis.set_label_coords(.04, .5)
		#rel_ax.set_xlim(left=1, right=10**9)
		rel_ax.set_ylabel("Relative Error", fontsize=10)
		rel_ax.legend(
			fontsize=6, 
			loc='upper left', 
			bbox_to_anchor=(.05, .95),
			labelspacing=1
		)
		rel_ax.tick_params(axis='both', which='major', labelsize=8)
		rel_ax.tick_params(axis='y', which='major', labelrotation=60)
		
		# Create offset transform by 5 points in x direction
		dx = 2/72.; dy = 9/72. 
		offset = matplotlib.transforms.ScaledTranslation(dx, dy, rel_fig.dpi_scale_trans)
		# apply offset transform to all x ticklabels.
		for label in rel_ax.yaxis.get_majorticklabels():
			label.set_transform(label.get_transform() + offset)

		#
		# ADDITIVE VERSION
		#
		
		# set labels and lims
		add_ax.set_xticks([10**x for x in range(10)])
		add_ax.set_xlabel("Rank", fontsize=10)
		add_ax.xaxis.set_label_coords(.5, .04)
		add_ax.yaxis.set_label_coords(.04, .5)
		#add_ax.set_xlim(left=1, right=10**9)
		add_ax.set_ylabel("Additive Error", fontsize=10)
		add_ax.legend(
			fontsize=6, 
			loc='upper left', 
			bbox_to_anchor=(.05, .95),
			labelspacing=1
		)
		add_ax.tick_params(axis='both', which='major', labelsize=8)
		add_ax.tick_params(axis='y', which='major', labelrotation=60)
		
		# Create offset transform by 5 points in x direction
		dx = 2/72.; dy = 9/72. 
		offset = matplotlib.transforms.ScaledTranslation(dx, dy, rel_fig.dpi_scale_trans)
		# apply offset transform to all x ticklabels.
		for label in add_ax.yaxis.get_majorticklabels():
			label.set_transform(label.get_transform() + offset)
		
		# save the graphs
		add_fig.savefig("add_log_temp.pdf")
		rel_fig.savefig("rel_log_temp.pdf")
		
		subprocess.run(["pdfcrop", "add_log_temp.pdf", "add_log.pdf"]) 
		subprocess.run(["pdfcrop", "rel_log_temp.pdf", "rel_log.pdf"])
		os.remove("add_log_temp.pdf")
		os.remove("rel_log_temp.pdf")
		
	def add_dataset(file, rel_ax, add_ax, colour):
		with open(file, mode='rb') as file:
			samples = pickle.load(file)
		points = samples.sample_points
		data = samples.data
		info = samples.info
		n = samples.n
		
		# plot additive
		x, y = Plotting.sparsify(points, [z/n for z in data.perc99])
		add, = add_ax.plot(x, y, c=colour, linestyle="-", lw=0.5)
		
		# plot relative
		x, y = Plotting.sparsify(points, [data.perc99[i]/max(points[i], 1) for i in range(len(points))])
		rel, = rel_ax.plot(x, y, c=colour, linestyle="-", lw=0.5)
		
		# set legend
		if "final" in info["user"]:
			legend = (
				f"Q: {info['Q']}\n"
				f"J: {info['J']}\n"
				f"capacity: {info['cap']}\n"
				f"improvement: {'YES' if info['improvement'] else 'NO'}"
			)
		else:
			legend = pformat(info)
		
		add.set_label(legend)
		rel.set_label(legend)
	
	def sparsify(x_in, y_in, s=50):
		x = []
		y = []
		i = 0
		while i <= len(x_in) - s:
			j, y_j = max(enumerate(y_in[ i : i+s ]), key=lambda t: t[1])
			x.append(x_in[i+j])
			y.append(y_j)
			i+=s
		return (x, y)

if __name__ == '__main__':
	main()
