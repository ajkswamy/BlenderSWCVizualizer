Ajayrama Kumaraswamy, 17.10.2014

BlenderSWCVizualizer
====================

Scripts for vizualizing a SWC files in Blender

'blenderHelper.py' contains classes for representing a SWC model in Blender

The other four files are examples of how to use this class. 

Mainly tested on a Ubuntu Linux Machine. Needs Blender to be installed.
 
Usages and Notes:
======

vizSWC.py: Used to visualize the model in the SWC file as it is.
----------
1. The argument triplet in "abc.importWholeSWC([0, 0, 1])" is the RGB value of the color with which the neuron will be colored
2. To use, change the string in sys.path.append('...') to the path which contains blenderHelper.py
3. Then, type in a terminal, "blender --python vizSWC.py -- /path/to/SWC/file".

compareNrnsBlender.py: Used to vizualize multiple SWC files together.
---
1. To use, specify the paths to the SWC files in the list 'swcs'.
2. Specify the colors with which the neurons are to be colored in 'cols'. The neuron to color mapping in element wise one to one.
3. Change the string in sys.path.append('...') to the path which contains blenderHelper.py
4. Then, in a terminal, enter "blender --python compareNrnsBlender.py"

vizSSWC.py: Used to vizualize the model in an SSWC file with possible coloring of segments.
-----
1. An SSWC file contains an extra column in addition to the 7 in a normal SWC file. This extra column must contain only integers from 0 to 7.
2. The array baseCols is a color map.
3. If a segment has the value 3 in it's 8th column, then it will be coloured with the 4th color in the colormap in baseCols.
4. To use, change the string in sys.path.append('...') to the path which contains blenderHelper.py
5. Then in a terminal, enter "blender --python vizSSWC.py -- /path/to/SWC/file"

vizSSWC_schematic.py: Similar to vizSSWC.py but with the radii of all the segments set to a single specified value.
----------
1. The radius is specified as the right hand side of "swcData[:, 5] = 0.01" 
2. To use, change the string in sys.path.append('...') to the path which contains blenderHelper.py
3. Then, in a terminal, enter "blender --python vizSSWC.py -- /path/to/SWC/file"

vizSSWCMultiple.py: Used to vizualize multiple SSWC files together.
---

1. To use, specify the paths to the SSWC files in the list 'swcFiles'
2. Change the string in sys.path.append('...') to the path which contains blenderHelper.py
3. Then, in a terminal, enter "blender --python vizSSWCMultiple.py"

Blender Basics:
===

1. Zoom View: Mouse wheel scroll
2. 3D rotate View: Mouse wheel click + drag
3. Select Object: Right click(hold Shift to select more)
4. Move Object: Select object, press G, move object and left click to fix position
5. Rotate Object: Select object, press R, rotate object and left click to fix position