# author: Stephen Lesage
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from scipy.linalg import norm
import matplotlib.pyplot as plt

# Setup figure
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection='3d')
ax.set_axis_off()
ax.set_xlim(-7, 7)
ax.set_ylim(-7, 7)
ax.set_zlim(-7, 7)

# Definitions
def xyz_to_r(r):
	# input: catesian x y z
	# output: spherical r
	return np.sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2])

# Can generate regular cones as well
def trunc_cone(length_i, dist_from_length_i, radius_i, radius_f, color, alpha, n_steps):
	# cone z_axis vector
	cone_zaxis = length_i
	# length of cone z_axis
	mag_zaxis = norm(cone_zaxis)
	# z-hat
	cone_zhat = cone_zaxis / mag_zaxis
	# find vector not along cone z-hat
	not_cone_zhat = np.array([1, 1, 0])
	if (cone_zhat == not_cone_zhat).all():
		not_cone_zhat = np.array([0, 1, 0])
	# make vector perp to cone z-hat
	perp_to_cone_zhat = np.cross(cone_zhat, not_cone_zhat)
	# x-hat
	cone_xhat = perp_to_cone_zhat / norm(perp_to_cone_zhat)
	# y-hat
	cone_yhat = np.cross(cone_zhat, cone_xhat)
	# dtheta-dz surface
	dlength = np.linspace(0, dist_from_length_i, n_steps)
	theta = np.linspace(0, 2 * np.pi, n_steps)
	# meshgrid
	dlength, theta = np.meshgrid(dlength, theta)
	# dradius
	dradius = np.linspace(radius_i, radius_f, n_steps)
	# create data
	x, y, z = [length_i[i] + cone_zhat[i] * dlength + dradius * np.sin(theta) * cone_xhat[i] + dradius * np.cos(theta) * cone_yhat[i] for i in [0, 1, 2]]
	# plot
	ax.plot_surface(x, y, z, color=color, alpha=alpha, linewidth=0, antialiased=False)

def trunc_cone_half(length_i, dist_from_length_i, radius_i, radius_f, color, alpha, n_steps):
	# cone z_axis vector
	cone_zaxis = length_i
	# length of cone z_axis
	mag_zaxis = norm(cone_zaxis)
	# z-hat
	cone_zhat = cone_zaxis / mag_zaxis
	# find vector not along cone z-hat
	not_cone_zhat = np.array([1, 1, 0])
	if (cone_zhat == not_cone_zhat).all():
		not_cone_zhat = np.array([0, 1, 0])
	# make vector perp to cone z-hat
	perp_to_cone_zhat = np.cross(cone_zhat, not_cone_zhat)
	# x-hat
	cone_xhat = perp_to_cone_zhat / norm(perp_to_cone_zhat)
	# y-hat
	cone_yhat = np.cross(cone_zhat, cone_xhat)
	# dtheta-dz surface
	dlength = np.linspace(0, dist_from_length_i, n_steps)
	theta = np.linspace(0, np.pi, n_steps)
	# meshgrid
	dlength, theta = np.meshgrid(dlength, theta)
	# dradius
	dradius = np.linspace(radius_i, radius_f, n_steps)
	# create data
	x, y, z = [length_i[i] + cone_zhat[i] * dlength + dradius * np.sin(theta) * cone_xhat[i] + dradius * np.cos(theta) * cone_yhat[i] for i in [0, 1, 2]]
	# plot
	ax.plot_surface(x, y, z, color=color, alpha=alpha, linewidth=0, antialiased=False)

def sphere(radius, color, alpha, n_steps):
	# dtheta dphi
	theta = np.linspace(0, 2 * np.pi, n_steps)
	phi = np.linspace(0, np.pi, n_steps)
	# manually generate caresian meshgrid
	x = radius * np.outer(np.cos(theta), np.sin(phi))
	y = radius * np.outer(np.sin(theta), np.sin(phi))
	z = radius * np.outer(np.ones(np.size(theta)), np.cos(phi))
	# plot
	ax.plot_surface(x, y, z, color=color, alpha=alpha, linewidth=0, antialiased=False)

def sphere_match_cone(pt_at_obj_origin, dist_to_obj_from_pt_on_plane, color, alpha, n_steps):
	# distance to object plane along zaxis
	dist_to_object_origin = norm(pt_at_obj_origin)
	# find sphere radius
	sphere_radius = np.sqrt(dist_to_obj_from_pt_on_plane * dist_to_obj_from_pt_on_plane + dist_to_object_origin * dist_to_object_origin)
	# dtheta tphi
	theta = np.linspace(0, 2 * np.pi, n_steps)
	phi = np.linspace(0, np.pi, n_steps)
	# manually generate caresian meshgrid
	x = sphere_radius * np.outer(np.cos(theta), np.sin(phi))
	y = sphere_radius * np.outer(np.sin(theta), np.sin(phi))
	z = sphere_radius * np.outer(np.ones(np.size(theta)), np.cos(phi))
	# plot
	ax.plot_surface(x, y, z, color=color, alpha=alpha, linewidth=0, antialiased=False)

def lens_in_cone(cone_i, cone_f, radius_i, radius_f, distance_away, color, alpha1, n_steps):
	# defining cone parameters
	dist_to_cone_beg = norm(cone_i)
	dist_to_cone_end = dist_to_cone_beg + cone_f
	cone_length = cone_f
	# defining lens prameters at cone edges
	dist_to_lens = dist_to_cone_beg + distance_away
	delta_r = radius_f - radius_i
	lens_radius = radius_i + (delta_r * float(distance_away) / float(cone_length))
	# find sphere radius
	sphere_radius = np.sqrt(lens_radius * lens_radius + dist_to_lens * dist_to_lens)
	# find angle of cone
	aperature = np.arctan(lens_radius / dist_to_lens)
	# shift spherical coordinates starting position
	new_phi_origin = ( np.arctan(np.sqrt(cone_i[0] * cone_i[0] + cone_i[1] * cone_i[1]) / cone_i[2]) )
	new_theta_origin = ( np.arctan(cone_i[1] / cone_i[0]) )
	# dtheta dphi
	theta = np.linspace(0, 2 * np.pi, n_steps)
	phi = np.linspace(0, aperature, n_steps)
	# manually generate caresian meshgrid
	x = np.zeros((n_steps,n_steps))
	y = np.zeros((n_steps,n_steps))
	z = np.zeros((n_steps,n_steps))
	# create data
	for i in range(0, n_steps):
		for j in range(0, n_steps):
			x[i][j] = sphere_radius * (-(np.sin(phi[i])*np.cos(new_phi_origin)*np.cos(new_theta_origin)*np.cos(theta[j])) - (np.sin(phi[i])*np.sin(new_theta_origin)*np.sin(theta[j])) + (np.cos(phi[i])*np.sin(new_phi_origin)*np.cos(new_theta_origin)))
			y[i][j] = sphere_radius * (-(np.sin(phi[i])*np.cos(new_phi_origin)*np.sin(new_theta_origin)*np.cos(theta[j])) + (np.sin(phi[i])*np.cos(new_theta_origin)*np.sin(theta[j])) + (np.cos(phi[i])*np.sin(new_phi_origin)*np.sin(new_theta_origin)))
			z[i][j] = sphere_radius * ( (np.sin(phi[i])*np.sin(new_phi_origin)                         *np.cos(theta[j]))                                                              + (np.cos(phi[i])*np.cos(new_phi_origin)                         ))
	# plot
	ax.plot_surface(x, y, z, color=color, alpha=alpha1, linewidth=0, antialiased=False)

def lens_dist_rad(lens_axis, lens_radius, color, alpha1, n_steps):
	# distance to lens origin
	dist_to_lens = norm(lens_axis)
	# find sphere radius
	sphere_radius = np.sqrt(lens_radius * lens_radius + dist_to_lens * dist_to_lens)
	# find angle of cone
	aperature = np.arctan(lens_radius / dist_to_lens)
	# shift spherical coordinates starting position
	new_phi_origin = ( np.arctan(np.sqrt(lens_axis[0] * lens_axis[0] + lens_axis[1] * lens_axis[1]) / lens_axis[2]) )
	new_theta_origin = ( np.arctan(lens_axis[1] / lens_axis[0]) )
	# dtheta dphi
	theta = np.linspace(0, 2 * np.pi, n_steps)
	phi = np.linspace(0, aperature, n_steps)
	# manually generate caresian meshgrid
	x = np.zeros((n_steps,n_steps))
	y = np.zeros((n_steps,n_steps))
	z = np.zeros((n_steps,n_steps))
	# create data
	for i in range(0, n_steps):
		for j in range(0, n_steps):
			x[i][j] = sphere_radius * ( (np.sin(phi[i])*np.cos(new_phi_origin)*np.cos(new_theta_origin)*np.cos(theta[j])) - (np.sin(phi[i])*np.sin(new_theta_origin)*np.sin(theta[j])) + (np.cos(phi[i])*np.sin(new_phi_origin)*np.cos(new_theta_origin)))
			y[i][j] = sphere_radius * (-(np.sin(phi[i])*np.cos(new_phi_origin)*np.sin(new_theta_origin)*np.cos(theta[j])) + (np.sin(phi[i])*np.cos(new_theta_origin)*np.sin(theta[j])) + (np.cos(phi[i])*np.sin(new_phi_origin)*np.sin(new_theta_origin)))
			z[i][j] = sphere_radius * ( (np.sin(phi[i])*np.sin(new_phi_origin)                         *np.cos(theta[j]))                                                              + (np.cos(phi[i])*np.cos(new_phi_origin)                         ))
	# plot
	ax.plot_surface(x, y, z, color=color, alpha=alpha1, linewidth=0, antialiased=False)

def lens_dist_angle(radius, angle, color, alpha, n_steps):
	# find sphere radius
	sphere_radius = xyz_to_r(radius)
	# angle of cone in rad
	aperature = angle * np.pi / 180
	# shift spherical coordinates starting position
	new_phi_origin = ( np.arctan(np.sqrt(lens_axis[0] * lens_axis[0] + lens_axis[1] * lens_axis[1]) / lens_axis[2]) )
	new_theta_origin = ( np.arctan(lens_axis[1] / lens_axis[0]) )
	# dtheta dphi
	theta = np.linspace(0, 2 * np.pi, n_steps)
	phi = np.linspace(0, aperature, n_steps)
	# manually generate caresian meshgrid
	x = np.zeros((n_steps,n_steps))
	y = np.zeros((n_steps,n_steps))
	z = np.zeros((n_steps,n_steps))
	# create data
	for i in range(0, n_steps):
		for j in range(0, n_steps):
			x[i][j] = sphere_radius * ( (np.sin(phi[i])*np.cos(new_phi_origin)*np.cos(new_theta_origin)*np.cos(theta[j])) - (np.sin(phi[i])*np.sin(new_theta_origin)*np.sin(theta[j])) + (np.cos(phi[i])*np.sin(new_phi_origin)*np.cos(new_theta_origin)))
			y[i][j] = sphere_radius * (-(np.sin(phi[i])*np.cos(new_phi_origin)*np.sin(new_theta_origin)*np.cos(theta[j])) + (np.sin(phi[i])*np.cos(new_theta_origin)*np.sin(theta[j])) + (np.cos(phi[i])*np.sin(new_phi_origin)*np.sin(new_theta_origin)))
			z[i][j] = sphere_radius * ( (np.sin(phi[i])*np.sin(new_phi_origin)                         *np.cos(theta[j]))                                                              + (np.cos(phi[i])*np.cos(new_phi_origin)                         ))
	# plot
	ax.plot_surface(x, y, z, color=color, alpha=alpha1, linewidth=0, antialiased=False)

##########################
########## MAIN ##########
##########################

pt1 = np.array([1, 1, 1])
dist_1 = 6
r1 = 1
r2 = 3
color = 'blue'
alpha = 0.5
n_steps = 30

#plot
#trunc_cone(pt1, dist_1, r1, r2, color, alpha, n_steps)
#trunc_cone(pt1, dist_1, r1, r2+1, color, alpha, n_steps)

trunc_cone(-pt1, dist_1, r1, r2, color, alpha, n_steps)
#trunc_cone(-pt1, dist_1, r1, r2+1, color, alpha, n_steps)

trunc_cone_half(pt1, dist_1, r1, r2, color, alpha, n_steps)

#sphere(xyz_to_r(pt1), color, alpha, n_steps)

sphere_match_cone(xyz_to_r(pt1), r1, color, alpha, n_steps)

#lens_dist_rad(pt1+2, r2, color, alpha, n_steps)

#lens_dist_angle(pt1, 30, color, alpha, n_steps)

lens_in_cone(pt1, dist_1, r1, r2, 1, color, alpha, n_steps)
lens_in_cone(pt1, dist_1, r1, r2, 2.5, color, alpha, n_steps)
lens_in_cone(pt1, dist_1, r1, r2, 3, color, alpha, n_steps)
lens_in_cone(pt1, dist_1, r1, r2, 4.5, color, alpha, n_steps)

plt.show()
