from ete2 import Tree

begoniaList = ['PLEisotig035802', 'CONisotig03447', 'CONisotig01392', 'PLEisotig03474', 'VENisotig01282', 'CONisotig017382', 'PLEisotig03031', 'VENisotig03204', 'CONisotig03330', 'PLEisotig00065', 'PLEisotig000712', 'PLEisotig000722', 'CONisotig03376']


f = open('RAxML_bestTree.begRosidsIpomoeaTree.raxml')
tree = f.read()
t = Tree(tree)


def printAllNodes(tree):
	#Gets all node names and returns a list of them
	nodeNames = []
	for node in tree.traverse('postorder'):
		nodeName = (node.name)
		if nodeName:
			nodeNames.append(nodeName)
	return(nodeNames)

def searchBySize(node, size):
	#Finds nodes with a given number of leaves
    matches = []
    for n in node.traverse():
       if len(n) == size:
       		matches.append(n)
    print(matches)


totalLeaves = []
for node in t.traverse('postorder'):
	leaves = node.iter_leaf_names()
	for i in leaves:
		totalLeaves.append(i)
totalLeaves = set(totalLeaves)

monophyleticList = []

for node in t.traverse('postorder'):	
	leafList = []
	leaves = node.iter_leaf_names()
	for l in leaves:
		leafList.append(l)
	
	if set(leafList) == set(begoniaList):
		print(leafList)



# this gives the node which contains all of the begonia sequences, in this case, the sequences in list begoniaList. When implementing in a pipeline get it to say whether there is one or not to output whether that gene family ismonophyletic for begonia





