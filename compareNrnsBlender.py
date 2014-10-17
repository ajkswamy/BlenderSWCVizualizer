#Ajayrama Kumaraswamy, 17.10.2014
import sys
import os
sys.path.append('/home/ajay/repos/BlenderSWCVizualizer')

from blenderHelper import BlenderSWCImporter


swcs = [

    'ExampleSWCs/HSN-fluoro01.CNG.swc',
    'ExampleSWCs/HSN-fluoro02.CNG.swc',
    'ExampleSWCs/HSN-fluoro03.CNG.swc',

       ]


# cols = [[ 0.        ,  0.        ,  0.5       ],
        # [ 0.        ,  0.00196078,  1.        ],
        # [ 0.        ,  0.50392157,  1.        ],
        # [ 0.08538899,  1.        ,  0.88235294],
        # [ 0.49019608,  1.        ,  0.47754586],
        # [ 0.89500316,  1.        ,  0.07273877],
        # [ 1.        ,  0.58169935,  0.        ],
        # [ 1.        ,  0.11692084,  0.        ]]

cols = [[1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
        [0.5, 0, 0.5],
        [0, 0.5, 0.5],
        [0.5, 0.5, 0]]

# cols = [[0, 0.2, 0.2], [0.1, 0, 0.1]]

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

    tmpB = BlenderSWCImporter(os.path.abspath(nrn), add, matchOrigin)
    tmpB.importWholeSWC(col=cols[nrnInd])

    nrnsBlender.append(tmpB)



