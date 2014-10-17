#Ajayrama Kumaraswamy, 17.10.2014
import sys
sys.path.append('/home/ajay/repos/BlenderSWCVizualizer')

from blenderHelper import BlenderSWCImporter
import numpy as np

print(sys.argv)
assert len(sys.argv) == 5, \
    'This script takes only 5 arguments, with the path of the swcfile expected as the 5th arguement, but ' \
    + str(len(sys.argv)) + ' found'
sswcFName = sys.argv[4]

swcData = np.loadtxt(sswcFName)
swcData[:, 5] = 0.01

nPts = 8

baseCols = np.asarray( [[ 0.        ,  0.        ,  0.5       ],
                        [ 0.        ,  0.00196078,  1.        ],
                        [ 0.        ,  0.50392157,  1.        ],
                        [ 0.08538899,  1.        ,  0.88235294],
                        [ 0.49019608,  1.        ,  0.47754586],
                        [ 0.89500316,  1.        ,  0.07273877],
                        [ 1.        ,  0.58169935,  0.        ],
                        [ 1.        ,  0.11692084,  0.        ]])

if nPts % len(baseCols) > 0:
    raise(ValueError('nPts must be a integral multiple of ' + str(len(baseCols))))

cols = np.zeros([nPts, 3])
scaleFactor = int(nPts /len(baseCols))

for ind in range(3):
    cols[:, ind] = np.interp(range(nPts), range(0, nPts, scaleFactor), baseCols[:, ind])



blenderObj = BlenderSWCImporter(sswcFName, colMap=cols, matchRootOrigin=False, swcData=swcData)
blenderObj.importWholeSWC()