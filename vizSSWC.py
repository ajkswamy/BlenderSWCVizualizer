import sys
sys.path.append('/home/aj/repos/BlenderSWCVizualizer')

from blenderHelper import BlenderSWCImporter
import numpy as np
import bpy

print(sys.argv)
assert len(sys.argv) == 5, \
    "Improper Usage! Please use as\n blender --python {} -- <path to SWC File>"
sswcFName = sys.argv[4]

# ----------------------------------------------------------------------------------------------------------------------
# --------------------------------Colors of SSWC Materials--------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------

# The range of values in the eight column of the input SSWC file will be divided into "nPts" number of bins.
# "cols", "alphas" and "emits" contain respectively the RGB triplet, transparency alpha and emit values of the material
# assoicated with these bins.
nPts = 8

baseCols = np.asarray( [[ 0.        ,  0.        ,  0.5       ],
                        [ 0.        ,  0.00196078,  1.        ],
                        [ 0.        ,  0.50392157,  1.        ],
                        [ 0.08538899,  1.        ,  0.88235294],
                        [ 0.49019608,  1.        ,  0.47754586],
                        [ 0.89500316,  1.        ,  0.07273877],
                        [ 1.        ,  0.58169935,  0.        ],
                        [ 1.        ,  0.11692084,  0.        ]])

# ----------------------------------------------------------------------------------------------------------------------
# ---------------------------Emit values of SSWC Materials--------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
emits = [1 for x in range(nPts)]
# ----------------------------------------------------------------------------------------------------------------------
# --------------------------Alpha values of SSWC Materials--------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
alphas = [1 for x in range(nPts)]
# alphas = [1, 1, 1, 0.5, 0.5, 1, 1, 1]
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
    mat.use_transparency = True
    mat.diffuse_color = scol
    mat.diffuse_intensity = 1.0
    mat.alpha = alphas[scolInd]
    mat.emit = emits[scolInd]
    sswcMaterials.append(mat)

blenderObj = BlenderSWCImporter(sswcFName, sswcMaterials=sswcMaterials, matchRootOrigin=False)
blenderObj.importWholeSWC()