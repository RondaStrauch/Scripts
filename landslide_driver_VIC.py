# ##!/usr/bin/env python
#############################################################################
"""
Purpose is to gather input data, run LandslideProability component, and
visualize/store data and outputs.  Input data is provide by the user and
consists of elevation from a DEM to provide topographic traits such as slope,
contributing area, and flow direction. User also supplies soil characteristics
derived from a soil survey, land cover, or other sources, including:
transmissivity, cohesion, internal friction angle, density, and soil thickness.

This demonstration uses groundwater recharge derived from the Variable
Infiltration Capacity (VIC) model and the 'data_driven_spatial' recharge
option of the component. A separate 'Source Tracking Algorithm' was run
to provide the draining 'Hydrology Source Domain' (HSD) grid IDs and fractions
used 'route' recharge to the model domain resolution of 30 m.

Output: Figures - relative wetness and probability of failure based on
Monte Carlo simulations of infinite slope factor-of-safety stability index.

See UPDATE for changing domains or hydrologic data model runs

@author: R.Strauch, E.Istanbulluoglu, S.S.Nudurupati, C.Bandaragoda
Univerity of Washington
Created on Thu Aug 20 16:47:11 2015
Last edit April 2017
"""
#############################################################################
# %% Import libraries and components
import numpy as np
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
# from landlab import RasterModelGrid
from landlab.plot.imshow import imshow_node_grid
from landlab.io import read_esri_ascii
from landlab.io import write_esri_ascii
from collections import defaultdict
from landlab.components.landslides import LandslideProbability
# from Re_Utilityv02 import calculate_recharge
import cPickle as pickle
# import pandas as pd

# %% load or create the data in fields and Set boundary conditions

# Location of parameter data
# UPDATE FOR DIFFERENT DOMAINS
data_folder = '/data1/landlab_eco/landslide_ronda/NOCA_analysis'

# load in DEM from ArcGIS
(grid, z) = read_esri_ascii(data_folder+'/elevation.txt',
                            name='topographic__elevation')
grid.at_node.keys()     # loads DEM grid with elevation
grid.set_nodata_nodes_to_closed(grid.at_node['topographic__elevation'], -9999)
# set boundary conditions closed where no data

# Load other fields from text files (ascii from ArcGIS)
# set_nodata_nodes_to_inactive that have no data, for each field

# load slope from text file
(grid1, slope) = read_esri_ascii(data_folder+'/slope_tang17.txt')
grid.add_field('node', 'topographic__slope', slope)
grid.set_nodata_nodes_to_closed(grid.at_node['topographic__slope'], -9999)
grid.set_nodata_nodes_to_closed(grid.at_node['topographic__slope'], 0.0)
# Load contributing area from text file
(grid1, ca) = read_esri_ascii(data_folder+'/cont_area.txt')
grid.add_field('node', 'topographic__specific_contributing_area', ca)
grid.set_nodata_nodes_to_closed(grid.at_node[
    'topographic__specific_contributing_area'], -9999)
# load transmissivity from text file
(grid1, T) = read_esri_ascii(data_folder+'/transmis_model.txt')
grid.add_field('node', 'soil__transmissivity', T)
grid.set_nodata_nodes_to_closed(grid.at_node['soil__transmissivity'], -9999)
# load cohesion from text file
(grid1, C) = read_esri_ascii(data_folder+'/cohesion_mean.txt')
C[C == 0.0] = 1.0  # ensure not 0 Pa for use in distributions generation
grid.add_field('node', 'soil__mode_total_cohesion', C)
grid.set_nodata_nodes_to_closed(grid.at_node[
    'soil__mode_total_cohesion'], -9999)
# (grid1, C_min) = read_esri_ascii(data_folder+'/cohesion_min.txt')
C_min = C*(0.3)
grid.add_field('node', 'soil__minimum_total_cohesion', C_min)
grid.set_nodata_nodes_to_closed(grid.at_node[
    'soil__minimum_total_cohesion'], -9999)
(grid1, C_max) = read_esri_ascii(data_folder+'/cohesion_max.txt')
grid.add_field('node', 'soil__maximum_total_cohesion', C_max)
grid.set_nodata_nodes_to_closed(grid.at_node[
    'soil__maximum_total_cohesion'], -9999)
# Load internal angle of friction from text file
(grid1, phi) = read_esri_ascii(data_folder+'/frict_angle.txt')
grid.add_field('node', 'soil__internal_friction_angle', phi)
grid.set_nodata_nodes_to_closed(grid.at_node[
    'soil__internal_friction_angle'], -9999)
# set soil density value and assign to all nodes
grid['node']['soil__density'] = 2000*np.ones(grid.number_of_nodes)
# load soil thickness from text file
(grid1, hs) = read_esri_ascii(data_folder+'/soil_depth_model.txt')
grid.add_field('node', 'soil__thickness', hs)
grid.set_nodata_nodes_to_closed(grid.at_node['soil__thickness'], -9999)

# load analysis mask and actual landslide from text file
(grid1, mask) = read_esri_ascii(data_folder+'/exclud_mask.txt')
grid.add_field('node', 'exclusion_mask', mask)
grid.set_nodata_nodes_to_closed(grid.at_node['exclusion_mask'], -9999)
(grid1, slides) = read_esri_ascii(data_folder+'/landslide_type.txt')
# 1-5 are landslides, 8 is no landslide mapped
grid.add_field('node', 'landslides', slides)
# load Frequency ratio probability from text file
(grid1, fr_prob) = read_esri_ascii(data_folder+'/fr_prob.txt')
grid.add_field('node', 'Fr__probability', fr_prob)

# %% number_of_iterations for Monte Carlo simulations
iterations = 100

# %% Recharge data Upload

# specify distribution to use in simulation
distribution = 'data_driven_spatial'
# Load pre-processed routed flows containing HSD_id and fractional drainage
# at each node and recharge dictionaries
    # !!! Ensure that the grid has been initialized with same
    # DEM (really) for pickle so same dimensions!!!
# dict of node id (key) and HSD_ids (values)
HSD_id_dict = pickle.load(open('dict_uniq_ids.p', 'rb'))
# dict of node id (key) and fractions (values)
fract_dict = pickle.load(open('dict_coeff.p', 'rb'))
# dict of HSD id (key) with arrays of recharge (values)
HSD_dict = pickle.load(open('HSD_dict.p', 'rb'))
#HSD ('hydrology source domain') dictionaries
HSD_inputs = [HSD_dict,HSD_id_dict, fract_dict]

# %% Istantiate and Run Infinite Slope model in Landslide Component
# Instantiate component for HSD distribution
LS_prob = LandslideProbability(grid,
    number_of_iterations=iterations,
    groudwater__recharge_distribution=distribution,
    groudwater__recharge_HSD_inputs=HSD_inputs)
# Run component
LS_prob.calculate_landslide_probability()

# %% Review outputs
#grid['node']['soil__mean_relative_wetness']
#grid['node']['landslide__mean_factor_of_safety']
#grid['node']['landslide__probability_of_failure']
#core_nodes = LS_prob.grid.core_nodes
##LS_prob.landslide__factor_of_safety_histogram[core_nodes[5036]]
# node_ids = range(0, grid.number_of_nodes)
# np.savetxt('node_ids.txt',node_ids)
# grid.at_node to see what fields are created
# range(0,grid.number_of_nodes) gives you a list of all node ids
# np.arange(0,grid.number_of_nodes).shape gives you how big this list is.
# grid.at_node['node_ids'] gives you back the field array

# %% Export data
# extract node id, landslide type, Fr probability, FS probability, slope
#sample_nodes = core_nodes
#data_extracted = {'prob_dyn': np.array(
#                 grid.at_node['landslide__probability_of_failure'][grid.core_nodes])}
#headers = ['prob_dyn']
#df = pd.DataFrame(data_extracted, index=sample_nodes, columns=(headers))
#df.to_csv('prob_dyn_n1000_modelsoil.csv')

# extract node ids for raster creation
#grid.add_field('node', 'node_ids', np.arange(0, grid.number_of_nodes))
#write_esri_ascii('node_id.txt',grid,names='node_ids')

# export probability of failure
# write_esri_ascii('prob_fail_lognorm_VIC.txt',grid,names='landslide__probability_of_failure')
    
# %% Figures & Plots

# plotting parameters
#mpl.rcParams['xtick.labelsize'] = 15
#mpl.rcParams['ytick.labelsize'] = 15
#mpl.rcParams['lines.linewidth'] = 1
#mpl.rcParams['axes.labelsize'] = 18
#mpl.rcParams['legend.fontsize'] = 15
#
#plt.clf
#plt.figure('Elevations from the DEM [m]')  # new fig, with a title
#imshow_node_grid(grid, 'topographic__elevation', cmap='terrain',
#                 grid_units=('coordinates', 'coordinates'),
#                 shrink=0.75, var_name='Elevation', var_units='m')
##elev = grid.node_vector_to_raster(z)
## control contour extent and placement of labels
##cs = pylab.contour(elev, extent=[-15,65535, -15,96495], hold='on', colors='gray')
##manual_locals = [(11000,11000), (14000,30000), (20000,8000), (36000,3000)]
##pylab.clabel(cs, inline=True, fmt='%1i', fontsize=10, manual=manual_locals)
#plt.savefig(outdir+'elevation.png', dpi=300)
#
#plt.figure('Mean Relative Wetness')
#imshow_node_grid(grid, 'soil__mean_relative_wetness', cmap='YlGnBu',
#                 grid_units=('coordinates', 'coordinates'),
#                 shrink=0.75, var_name='Relative Wetness',
#                 var_units='no units')
##cs = pylab.contour(elev, extent=[0,42500, 0,32500], hold='on', colors='gray')
##pylab.clabel(cs, inline=True, fmt='%1i', fontsize=10, manual=manual_locals)
#plt.savefig(outdir+'wetness.png', dpi=300)
#
#plt.figure('Probability of Failure')
#imshow_node_grid(grid, 'landslide__probability_of_failure', cmap='OrRd',
#                 grid_units=('coordinates', 'coordinates'), shrink=0.75,
#                 var_name='Probability of Failure', var_units='no units')
##cs = pylab.contour(elev, extent=[0,42500, 0,32500], hold='on', colors='black')
##pylab.clabel(cs, inline=True, fmt='%1i', fontsize=10, manual=manual_locals)
#plt.savefig(outdir+'Prob_failure.png', dpi=300)

#example_FS_dist = LS_prob.landslide__factor_of_safety_histogram[
#    core_nodes[48]]
#plt.figure('Ex_FS_Distribution')
#plt.hist(example_FS_dist, 22)
#plt.plot([1, 1], [0, 10], color='r', linestyle='-', linewidth=3)
#plt.title('Sample FS Distribution at one node')
#plt.ylabel('Frequency')
#plt.xlabel('Calculated Factor-of-Safety')
#plt.savefig(outdir+'ex_FShist.png', dpi=300)
#
#grid['node']['landslide__probability_of_failure'][core_nodes[48]]
