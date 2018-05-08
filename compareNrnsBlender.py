import sys
import os
import numpy as np
sys.path.append('/home/aj/repos/BlenderSWCVizualizer')

from blenderHelper import BlenderSWCImporter
import bpy

swcs = [

    'ExampleSWCs/HSN-fluoro01.CNG.swc',
    'ExampleSWCs/HSN-fluoro02.CNG.swc',
    'ExampleSWCs/HSN-fluoro03.CNG.swc',

       ]

# ----------------------------------------------------------------------------------------------------------------------
# --------------------------------Colors of SWC Materials--------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------

cols = [[ 0.        ,  0.        ,  0.5       ],
        [ 0.        ,  0.00196078,  1.        ],
        [ 0.        ,  0.50392157,  1.        ],
        [ 0.08538899,  1.        ,  0.88235294],
        [ 0.49019608,  1.        ,  0.47754586],
        [ 0.89500316,  1.        ,  0.07273877],
        [ 1.        ,  0.58169935,  0.        ],
        [ 1.        ,  0.11692084,  0.        ]]


# ----------------------------------------------------------------------------------------------------------------------
# --------------------------------Colors of SSWC Materials--------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------

nPts = 8


SbaseCols = np.array([
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
# ---------------------------Emit values of SSWC Materials--------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
emits = [1 for x in range(nPts)]
# ----------------------------------------------------------------------------------------------------------------------
# --------------------------Alpha values of SSWC Materials--------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
alphas = [1 for x in range(nPts)]
# alphas = [1, 1, 1, 0.5, 0.5, 1, 1, 1]
# ----------------------------------------------------------------------------------------------------------------------

if nPts % len(SbaseCols) > 0:
    raise(ValueError('nPts must be a integral multiple of ' + str(len(SbaseCols))))

Scols = np.zeros([nPts, 3])
scaleFactor = int(nPts / len(SbaseCols))

for ind in range(3):
    Scols[:, ind] = np.interp(range(nPts), range(0, nPts, scaleFactor), SbaseCols[:, ind])

# ----------------------------------------------------------------------------------------------------------------------
nrnsBlender = []
matchOrigin = False


for nrnInd, nrn in enumerate(swcs):

    if nrnInd == 0:
        add = False
    else:
        add = True

    # tmpB = BlenderSWCImporter(nrn, add, matchOrigin, colMap=cols)
    # tmpB.importWholeSWC()

    # matchOrigin = True

    sswcMaterials = []
    for scolInd, scol in enumerate(Scols):
        mat = bpy.data.materials.new("Material {}".format(scolInd))
        mat.use_transparency = True
        mat.diffuse_color = scol
        mat.diffuse_intensity = 1.0
        mat.alpha = alphas[scolInd]
        mat.emit = emits[scolInd]
        sswcMaterials.append(mat)

    tmpB = BlenderSWCImporter(os.path.abspath(nrn), add, matchOrigin, sswcMaterials=sswcMaterials)
    tmpB.importWholeSWC(col=cols[nrnInd])

    nrnsBlender.append(tmpB)




