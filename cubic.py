import logging	
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle, PathPatch
import mpl_toolkits.mplot3d.art3d as art3d
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from itertools import product, combinations, combinations
from random import uniform, randint
from scipy.integrate import dblquad
from scipy.spatial import ConvexHull, Delaunay

class Cuboid(object):
	def __init__(self, o, length_=1, width_=1, height_=1):
		
		self.dimension = 3 # TODO n
		self.origin = o
		self.length = length_
		self.width = width_
		self.height = height_

	def print_cuboid(self):
		print (' origin = ', self.origin, '\n')
		print (' length = ', self.length, '\n' )
		print (' width  = ', self.width , '\n')
		print (' height = ', self.height, '\n')

	def draw(self, ax, o, color = 'b', alph = 0.1, centers=False):
		# Draw the cube based on a point (x, y, z), height, length and width.			

		x, y, z_ = o[0], o[1], o[2]

		if centers:
			ax.plot([x], [y], [z_], 'ro')
			ax.plot([x+self.length/2], [y+self.width/2], [z_+self.height/2], 'ro')
			ax.plot([x+self.length/2,x+self.length/2], [y+self.width/2,y+self.width/2], [z_, self.height], 'r-')
		
		if False:
			ax.plot([x, x], [y,y], [0,1], 'k-')
			ax.plot([x, x], [y+self.height,y+self.height], [0, self.height], 'm-')
			ax.plot([x + self.length,x+self.length], [y,y], [0, self.height], 'r-')
			ax.plot([x + self.length, x+self.length], [y+self.height, y+self.height], [0, self.height], 'g-')
	
		height_, length_ = self.height, self.width
		side = Rectangle((y, z_), length_, height_, facecolor=color, alpha=alph)
		ax.add_patch(side)
		art3d.pathpatch_2d_to_3d(side, z=x, zdir='x')
	
		height_, length_ = self.height, self.length 	
		side = Rectangle((x, z_), length_, height_, facecolor=color, alpha=alph)
		ax.add_patch(side)
		art3d.pathpatch_2d_to_3d(side, z=y, zdir='y')
	
		height_, length_ = self.width, self.length
		side = Rectangle((x, y), length_, height_, facecolor=color, alpha=alph)
		ax.add_patch(side)
		art3d.pathpatch_2d_to_3d(side, z=z_, zdir='z')
	
		height_, length_ = self.height, self.width 
		side = Rectangle((y, z_), length_, height_, facecolor=color, alpha=alph)
		ax.add_patch(side)
		art3d.pathpatch_2d_to_3d(side, z=x+self.length, zdir='x')
	
		height_, length_ = self.height, self.length 
		side = Rectangle((x, z_), length_, height_, facecolor=color, alpha=alph)
		ax.add_patch(side)
		art3d.pathpatch_2d_to_3d(side, z=self.width+y, zdir='y')
	
		height_, length_ = self.width, self.length
		side = Rectangle((x, y), length_, height_, facecolor=color, alpha=alph)
		ax.add_patch(side)
		art3d.pathpatch_2d_to_3d(side, z=z_+self.height, zdir='z')

#======================================================================================================================================================================================
# The function
def f(x, y, rho_x, rho_y, mu_x, mu_y, zeta_x, zeta_y, beta, delta_x, delta_y, gamma):
	exponent = delta_x - np.power((zeta_x*x - mu_x), rho_x) + delta_y - np.power((zeta_y*y - mu_y), rho_y)
	return gamma + beta * np.exp(exponent)


#======================================================================================================================================================================================
# Fitting the function to the cube

def funcuboid(cuboid, rho, ax, plot=False, view_squares=False, find_volume=False, function_alpha=0.1):
	
	# init
	ox, oy, oz = [cuboid.origin[_] for _ in range(cuboid.dimension)]		
	length, width, height = cuboid.length, cuboid.width, cuboid.height
	o = [ox, oy, oz]
	
	# filling

	delta_x,  delta_y  = 0, 0  # unused
	delta = [delta_x,  delta_y]
	beta = height
	gamma = oz  # to be initialized to the coords of the rectangle base (z)

	zeta_x  = 2. / ((rho/(rho-1))*length)	
	zeta_y  = 2. / ((rho/(rho-1))*width)	
	zeta = [zeta_x, zeta_y]
			
	# cube center
	cx = ox + length/2
	cy = oy + width/2
	
	mu_x = cx * zeta_x
	mu_y = cy * zeta_y
	mu = [mu_x, mu_y]

	if False:
		ax.plot([1, 1], [-2, -2], [0,1], 'ro-')

	theta_x =  (1/zeta_x) * np.power( (rho-1)/rho, 1/rho )	# Lx/2
	theta_y =  (1/zeta_y) * np.power( (rho-1)/rho, 1/rho )  # Ly/2
	
	x_3 =  mu_x/zeta_x - theta_x
	x_1 =  mu_x/zeta_x + theta_x
	
	x_4 =  mu_y/zeta_y - theta_y	
	x_2 =  mu_y/zeta_y + theta_y		
	
	corners = [[x_3, x_1], [x_4, x_2]]
	
	if plot and view_squares:
		ax.plot([x_1, x_1], [x_2, x_2], [0,beta], 'b-')
		ax.plot([x_1, x_1], [x_4, x_4], [0,beta], 'b-')
		ax.plot([x_3, x_3], [x_4, x_4], [0,beta], 'b-')
		ax.plot([x_3, x_3], [x_2, x_2], [0,beta], 'b-')

		ax.plot([x_1, x_1, x_3, x_3, x_1 ], [x_2, x_4, x_4, x_2, x_2], [0,0,0,0,0], 'b+-')
		ax.plot([x_1, x_1, x_3, x_3, x_1 ], [x_2, x_4, x_4, x_2, x_2], [beta, beta, beta, beta, beta], 'b+-')


		ax.plot([x_1], [x_2], 'm^', label='2_x') # A
		ax.plot([x_1], [x_4], 'b^') # B
		ax.plot([x_3], [x_4], 'b^') # C
		ax.plot([x_3], [x_2], 'b^', label='1_x') # D
				
	if plot:
		# central axis of the function curve
		ax.plot([mu_x/zeta_x, mu_x/zeta_x], [mu_y/zeta_y, mu_y/zeta_y], [0, beta], 'r-', label='$\mu/\zeta$')

		# start point of cube	
		ax.plot([ox], [oy], [0], 'r^', label='$r_x, r_y$')

		# start point of function curve
		ax.plot([x_3], [x_4], [0], 'ro', label='$s$')

	return o, beta, delta, gamma, zeta, mu, corners

#======================================================================================================================================================================================
def plot_function(ax, rho, mu, zeta, beta, delta, gamma, minmax, function_alpha, linewidth_):				
	min_x, max_x, step = minmax[0], minmax[1], minmax[2]
	x = y = np.arange(min_x, max_x, step)
	X, Y = np.meshgrid(x, y)
	zs = np.array([f(x, y, rho, rho, mu[0], mu[1], zeta[0], zeta[1], beta, delta[0], delta[1], gamma) for x,y in zip(np.ravel(X), np.ravel(Y))])
	Z = zs.reshape(X.shape)			
	ax.plot_surface(X, Y, Z, alpha=function_alpha, linewidth=linewidth_, edgecolors='b')

#======================================================================================================================================================================================
def paramprint(o, beta, delta, rho, gamma, zeta, mu, corners):	
	print ('\n       o = ', o)
	print ('\n    beta = ', beta)
	print ('\n   delta =', delta)
	print ('\n   rho   =', rho)
	print ('\n   gamma = ', gamma)
	print ('\n    zeta = ', zeta)
	print ('\n      mu = ', mu)
	print ('\n corners = ', corners)
	print ('\n________________________________________\n')
#======================================================================================================================================================================================

def main(rho, function_alpha, cube_alpha, view_cube, view_squares, find_volume, folder='figures/'):	
	fig = plt.figure(figsize = (14, 12))
	ax = fig.gca(projection='3d')
	ax.set_aspect("auto")
	ax.set_autoscale_on(True)
	
	_lim = 10.

	logger = logging.getLogger('myapp')
	hdlr = logging.FileHandler('./results1.csv')
	formatter = logging.Formatter('%(message)s')
	hdlr.setFormatter(formatter)
	logger.addHandler(hdlr) 
	logger.setLevel(logging.WARNING)
	logger.error('PointsInHull, PointsInConcavity, VInteg, VCube')	
		
	# draw the cube(s)
	nb_cubes = 1
	for cube in range(nb_cubes):	
		print (' ############')
		print ('  Cube #', cube)
		print (' ############')
		ox, oy, oz = uniform(-10, 10), uniform(-10, 10), 0
		o = (ox, oy, oz)
		length, width, height = uniform(1, 4), uniform(1, 6), uniform(1, 3)
		
		C = Cuboid(o, length_=length, width_=width, height_=height)
		C.draw(ax, o, color = 'r', alph = 0.1, centers=False)
		
		o, beta, delta, gamma, zeta, mu, corners = funcuboid(C, rho, ax, plot=True, view_squares=False, find_volume=False, function_alpha=0.01)
		paramprint(o, beta, delta, rho, gamma, zeta, mu, corners)	
		plot_function(ax, rho, mu, zeta, beta, delta, gamma, [-_lim, _lim, .1], function_alpha=0.15, linewidth_=0.1)			

		# volume part
		'''
		Checking if randommly generated points fall in both spaces.
		Assumptions of concavity of the constraints, which allows us to check for appartenance.
		
		using a function form for the constraints has the advantage of simplyfing the computation of the utility of a contract,
		by returing its weights (embed it in the function f, i.e. using the beta? gamma? to represent the w_k, weight of constraint c_k).
		the sutility will be the sum over the functions (constraints.)		
		'''
		#{{
		
		x_3, x_1 = corners[0]
		x_4, x_2 = corners[1]
		height = beta
		mu_x, mu_y = mu
		zeta_x, zeta_y = zeta
		delta_x, delta_y = delta

		n_points =  randint(20, 100)
		uniform_points = False
		
		points = [[0,0,0]] * n_points
		for i in range(n_points):
			if uniform_points:	# 1. Uniform, points anywhere
				points[i] = [uniform(min_x, max_x), uniform(min_x, max_x), uniform(0, zl)]
			
			else:	# 2. In the cube
				points[i] = [uniform(x_3, x_1), uniform(x_4, x_2), uniform(0, height)]
			# the contract point
			ax.scatter([points[i][0]], [points[i][1]], [points[i][2]], marker='o', c='c', s=1)
	
		points = np.array(points)
		
		# Construct the convex hull of the cube
		hull = ConvexHull(points)
		
		def in_hull(P, H):
			# Test if points in P are in H. P should be a n*k coordinates of n points in k dimension
			# H is either a scipy.spatial.Delaunay object or the m*k array of the coordinates of m points in k-dimension for which a Delaunay triangulation will be computed
			if not isinstance(H, Delaunay):
				H = Delaunay(H)
			return H.find_simplex(P)>=0	
		
		delaunay = Delaunay(points)
		n_points_in_cubes_hull = 0
		n_points_in_fs_concavity = 0
		for i in range(n_points):
			# check if p[i] is in the cube's hull
			if in_hull(points[i], delaunay):
				n_points_in_cubes_hull += 1
	
			# check if p[i] is in the function concavity
			f_of_p = f(points[i][0], points[i][1], rho, rho, mu_x, mu_y, zeta_x, zeta_y, beta, delta_x, delta_y, gamma)
			if f_of_p > 0:
				n_points_in_fs_concavity += 1
		
		print ('\nPoints in cube hull : ', n_points_in_cubes_hull, '/', n_points)
		print ('          concavity : ', n_points_in_fs_concavity, '/', n_points)
		#}}
		
		print ('\nComparing Volumes:')
			
		def integrand(y, x):
		    'y must be the first argument, and x the second.'
		    return f(x, y, rho, rho, mu_x, mu_y, zeta_x, zeta_y, beta, delta_x, delta_y, gamma)
		vinteg, err = dblquad(integrand,	# http://kitchingroup.cheme.cmu.edu/blog/2013/02/02/Integrating-functions-in-python/
				   x_3, x_1,		
				   lambda x: x_4, lambda x: x_2)
			
		print ('\t V_integ = ', vinteg)
		vcube = height * length * width
		print ('\t V_cube  = ', vcube)
	
		# PointsInHull, PointsInConcavity, VInteg, VCube
		logger.error('%f, %f, %f, %f' % ( n_points_in_cubes_hull, n_points_in_fs_concavity, vinteg, vcube ))	

	# TODO compute the contracts utlity using the cube (utility = beta if x in cube) and using f and show the equivalence.
	
	ax.set_xlabel('x (length)')
	ax.set_ylabel('y (width)')
	ax.set_zlabel('z=f(x,y), height')
	plt.xticks(np.arange(-_lim, _lim, 2))
	plt.yticks(np.arange(-_lim, _lim, 2))
	#zl = gamma + beta * 2.
	#ax.set_zlim([0., zl])
	ax.set_zlim([0., 3.1])
	plt.legend(fontsize=16)
	plt.title(r'$f(x, y; \rho, \beta, \gamma, \mu, \zeta)= \gamma+\beta e^{-(\zeta_1 x-\mu_1)^\rho-(\zeta_2 y-\mu_2)^\rho}$', fontsize=26)	
	
	plt.savefig(folder + 'Cube.pdf', format='pdf', dpi=1000)
	


	
if __name__ == '__main__':
	
	
	main(5000,
	     function_alpha = .24,
	     cube_alpha     = .04,
	     view_cube      = True,
	     view_squares   = False,
	     find_volume    = False)



