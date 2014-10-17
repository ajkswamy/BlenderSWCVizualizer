#Ajayrama Kumaraswamy, 17.10.2014
import sys
sys.path.append('/home/ajay/repos/BlenderSWCVizualizer')

from blenderHelper import BlenderSWCImporter

print(sys.argv)
assert len(sys.argv) == 5, \
    'This script takes only 5 arguments, with the path of the swcfile expected as the 5th arguement, but ' \
        + str(len(sys.argv)) + ' found'
swcFName = sys.argv[4]

abc = BlenderSWCImporter(swcFName, matchRootOrigin=False)
abc.importWholeSWC([0, 0, 1])