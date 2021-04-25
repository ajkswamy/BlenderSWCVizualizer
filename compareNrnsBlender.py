import sys
import os
import numpy as np
sys.path.append('/home/ajay/repos/BlenderSWCVizualizer')

from blenderHelper import BlenderSWCImporter
import bpy

import pathlib as pl

swcs = [

    'ExampleSWCs/HSN-fluoro01.CNG.swc',
    'ExampleSWCs/HSN-fluoro02.CNG.swc',
    'ExampleSWCs/HSN-fluoro03.CNG.swc',
       ]

# ----------------------------------------------------------------------------------------------------------------------
# --------------------------------Colors of SWC Materials--------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------

# # tab20 colormap
cols = [
    [0.12156863, 0.46666667, 0.70588235, 1.        ],
    [0.68235294, 0.78039216, 0.90980392, 1.        ],
    [1.        , 0.49803922, 0.05490196, 1.        ],
    [1.        , 0.73333333, 0.47058824, 1.        ],
    [0.17254902, 0.62745098, 0.17254902, 1.        ],
    [0.59607843, 0.8745098 , 0.54117647, 1.        ],
    [0.83921569, 0.15294118, 0.15686275, 1.        ],
    [1.        , 0.59607843, 0.58823529, 1.        ],
    [0.58039216, 0.40392157, 0.74117647, 1.        ],
    [0.77254902, 0.69019608, 0.83529412, 1.        ],
    [0.54901961, 0.3372549 , 0.29411765, 1.        ],
    [0.76862745, 0.61176471, 0.58039216, 1.        ],
    [0.89019608, 0.46666667, 0.76078431, 1.        ],
    [0.96862745, 0.71372549, 0.82352941, 1.        ],
    [0.49803922, 0.49803922, 0.49803922, 1.        ],
    [0.78039216, 0.78039216, 0.78039216, 1.        ],
    [0.7372549 , 0.74117647, 0.13333333, 1.        ],
    [0.85882353, 0.85882353, 0.55294118, 1.        ],
    [0.09019608, 0.74509804, 0.81176471, 1.        ],
    [0.61960784, 0.85490196, 0.89803922, 1.        ]
]


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

if nPts % len(SbaseCols) > 0:
    raise(ValueError('nPts must be a integral multiple of ' + str(len(SbaseCols))))

Scols = np.zeros([nPts, 4])
scaleFactor = int(nPts / len(SbaseCols))

for ind in range(4):
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
        mat.diffuse_color = scol
        sswcMaterials.append(mat)

    tmpB = BlenderSWCImporter(os.path.abspath(nrn), add, matchOrigin, sswcMaterials=sswcMaterials)
    tmpB.importWholeSWC(col=cols[nrnInd])

    nrnsBlender.append(tmpB)




