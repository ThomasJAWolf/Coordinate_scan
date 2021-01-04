"""
Code example for coordinate scans along the internal coordinates of a molecule. In this example, the 
photochemical ring-opening motion of alpha-terpinene is used.
Created by Thomas Wolf, 01/03/2021
"""

import numpy as np
from Diffraction import mol_geom # available from https://github.com/ThomasJAWolf/Diffraction_simulation
from Coordinate_scan import Coordinate_scan
   
############################################################################################################
## Example code ############################################################################################
############################################################################################################
# Numbering scheme (ordering of lines in start.xyz) for start geometry (ring-closed alpha-terpinene):
#
#         H18 H19 H20
#            \ | /
#             C8
#             |
#        H14  C5 H9
#          \ / \/
#          C6  C4-H10
#          |   |
#         C1  C3-H11
#        / \ / \
#     H13   C2 H12
#           |
#          C7
#        / | \
# H26-C17 H15 C16-H22
#    / \     / \
#  H25 H24  H23 H21

# Numbering scheme (ordering of lines in end.xyz) for end geometry (ring-opened alpha-terpinene):
#
#         H19 H18 H17
#            \ | /
#             C8
#             |
#        H16  C7 H13
#          \ / \/
#          C6  C9-H14
#          |   
#         C5  C1-H12
#        / \ / \
#     H15   C2 H11
#           |
#          C3
#        / | \
# H23-C10  H20 C4-H24
#    / \     / \
#  H22 H21  H26 H25

# Numbering scheme after reordering. It is important to anticipate the open ring connectivity in the numbering
# scheme of the start geometry. Otherwise, the interpolation will not yield chemically meaningful results.
#
#         H17 H16 H15
#            \ | /
#             C7
#             |
#        H14  C5 H11
#          \ / \/
#          C4  C6-H12
#          |   |
#         C3  C1-H10
#        / \ / \
#     H13   C2 H9
#           |
#          C8
#        / | \
# H26-C20 H18 C19-H22
#    / \     / \
#  H25 H24  H23 H21

# Reorder geometry according to numbering scheme (geomreord and elemreord) and define z-matrix connectivity

# Mapping of the order in start.xyz on the reordered geometry (numbers from 0 to number of atoms minus one).
Startnumbering = np.array([3,2,1,6,5,4,8,7,12,11,9,10,13,14,18,20,19,15,17,16,24,26,25,22,21,23])-1
# Connectivity (z-matrix style) of the reordered geometry.
Connectivity = [[],[0],[1,0],[2,1,0],[3,2,1],[4,3,2],[4,3,2],[1,2,0],[0,1,2],[0,1,2],[5,4,3],[5,4,3],[2,1,0],[3,2,1], \
                [6,1,0],[6,1,14],[6,1,14],[7,1,0],[7,1,17],[7,1,17],[18,7,1],[18,7,20],[18,7,20],[19,7,1],[19,7,23], \
                [19,7,23]]
# Mapping of the order in end.xyz on the reordered geometry (numbers from 0 to number of atoms minus one).
Endnumbering = np.array([1,2,5,6,7,9,8,3,11,12,13,14,15,16,17,18,19,20,4,10,24,25,26,21,22,23])-1
# Start geometry:
startgeo = mol_geom('start.xyz')
# End geometry:
endgeo = mol_geom('end.xyz')
terp = Coordinate_scan(startgeo,endgeo,Startnumbering, \
                          Endnumbering,Connectivity)
# Interpolate 181 points (including start and end geometry) between the start and end geometry.
terp.interpolate(181)
# Write interpolation result as xyz file:
terp.writexyzs(np.arange(181),'test')
