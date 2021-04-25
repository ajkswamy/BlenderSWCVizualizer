import sys
sys.path.append('/home/aj/repos/BlenderSWCVizualizer')

from blenderHelper import BlenderSWCImporter
import numpy as np
import bpy

print(sys.argv)
assert len(sys.argv) == 5, \
    "Improper Usage! Please use as\n blender --python {} -- <path to SWC File>"
sswcFName = sys.argv[4]
swcData = np.loadtxt(sswcFName)
swcData[:, 5] = 0.01


# ----------------------------------------------------------------------------------------------------------------------
# --------------------------------Colors of SSWC Materials--------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------

# The range of values in the eight column of the input SSWC file will be divided into "nPts" number of bins.
# "cols" contains the RGBA valuesof the material assoicated with these bins.
nPts = 8

baseCols = np.array([
    [ 0.78235294,  0.        ,  0.        ,  1.        ],
    [ 1.        ,  0.14509804,  0.14509804,  1.        ],
    [ 1.        ,  0.70980392,  0.70980392,  1.        ],
    [ 0.73333333,  0.89019608,  0.03921569,  0.5       ],
    [ 0.73333333,  0.89019608,  0.03921569,  0.5       ],
    [ 0.70980392,  0.70980392,  1.        ,  1.        ],
    [ 0.14509804,  0.14509804,  1.        ,  1.        ],
    [ 0.        ,  0.        ,  0.69529412,  1.        ]
    ]
    )

# ----------------------------------------------------------------------------------------------------------------------

if nPts % len(baseCols) > 0:
    raise(ValueError('nPts must be a integral multiple of ' + str(len(baseCols))))

cols = np.zeros([nPts, 3])
scaleFactor = int(nPts / len(baseCols))

for ind in range(3):
    cols[:, ind] = np.interp(range(nPts), range(0, nPts, scaleFactor), baseCols[:, ind])

sswcMaterials = []
for scolInd, scol in enumerate(cols):
    mat = bpy.data.materials.new("Material {}".format(scolInd))
    mat.diffuse_color = scol
    sswcMaterials.append(mat)

blenderObj = BlenderSWCImporter(sswcFName, sswcMaterials=sswcMaterials, matchRootOrigin=False, swcData=swcData)
blenderObj.importWholeSWC()