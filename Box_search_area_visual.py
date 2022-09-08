'''
Modified by Josh Tseng 2/28/2022
This is a script which can help you visuallize the search area of the structrues for AutoDock vina. You will have to run this script in pymol. And click on the atom of you wish to be the center atom. By executing the center_atom_coord command, it will print out the coordinates of the atom selected. By executing the Box command, it will draw a psedoatom and the area specified. 
Examples:
1. Get the coordinates of the center atom.
center_atom_coord('sele')

2. Box the search area. 
Box('sele',6,6,6)
'''


def center_atom_coord(selection):
	'''
	By default it uses the CA atom of the ligand as the central atom. This can be varied according to different ligand
	'''
	cen_coords = cmd.get_coords(selection,1)[0]
	print('Selected central coordinates is:', cen_coords)

	return cen_coords

def Box(selection, x, y, z):
	'''
	Run the script in pymol to visualize the search area.
	Sets up the search box within the protein, which is
	used in the docking protocol.
	'''
	
	pX, pY, pZ = center_atom_coord(selection)
	pymol.cmd.pseudoatom('Position', pos=[pX, pY, pZ])
	([X, Y, Z],[a, b, c]) = pymol.cmd.get_extent('Position')
	pymol.cmd.show('spheres', 'Position')
	minX = X+float(x)
	minY = Y+float(y)
	minZ = Z+float(z)
	maxX = X-float(x)
	maxY = Y-float(y)
	maxZ = Z-float(z)
	boundingBox = [BEGIN, LINES,
		VERTEX, minX, minY, minZ,
		VERTEX, minX, minY, maxZ,
		VERTEX, minX, maxY, minZ,
		VERTEX, minX, maxY, maxZ,
		VERTEX, maxX, minY, minZ,
		VERTEX, maxX, minY, maxZ,
		VERTEX, maxX, maxY, minZ,
		VERTEX, maxX, maxY, maxZ,
		VERTEX, minX, minY, minZ,
		VERTEX, maxX, minY, minZ,
		VERTEX, minX, maxY, minZ,
		VERTEX, maxX, maxY, minZ,
		VERTEX, minX, maxY, maxZ,
		VERTEX, maxX, maxY, maxZ,
		VERTEX, minX, minY, maxZ,
		VERTEX, maxX, minY, maxZ,
		VERTEX, minX, minY, minZ,
		VERTEX, minX, maxY, minZ,
		VERTEX, maxX, minY, minZ,
		VERTEX, maxX, maxY, minZ,
		VERTEX, minX, minY, maxZ,
		VERTEX, minX, maxY, maxZ,
		VERTEX, maxX, minY, maxZ,
		VERTEX, maxX, maxY, maxZ,
		END]
	boxName = 'Box'
	pymol.cmd.load_cgo(boundingBox, boxName)
	return(boxName)



