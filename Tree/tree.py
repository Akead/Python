#!/usr/bin/env python3

def author():
	#Paweł Lipiór
	print('\n' + 30*'%' + '\nAuthor: Paweł Lipiór \n' + 30*'%')

#------------------------------------- NODE --------------------------------------------------------------------

class Node:
	#Class representing one node of binary tree
	
	def __init__(self, key = None, left = None, right = None): 
		self.key = key
		self.left = left
		self. right = right

	def __str__(self):
		return str(self.key)

#------------------------------------- NODE PROP -----------------------------------------------------------------

class Node_Prop:
	#Class to perform basic operation on Node
	
	def __init__(self, node):
		self.node = node

	def sum_of_binary_tree(self):
		#Function to recursively sum nodes keys of binary tree
		def h_function(node):
				if (node == None):
					return 0
				return (node.key  + h_function(node.left) + h_function(node.right))
		return h_function(self.node)
		

	def avg_of_binary_tree(self):
		#Function to calculate averange value of nodes in tree
		def h_function(node, sum_of_tree = 0, count = 0):
			if node is None:
				return (0 ,0)
			else:
				(l_sum, l_count) = h_function(node.left)
				(r_sum, r_count) = h_function(node.right)
				return (node.key + l_sum + r_sum, 1 + l_count + r_count)
		(sum_of_tree, count) = h_function(self.node)
		return sum_of_tree/count if count > 0 else None


	def median_of_binary_tree(self):
		#Function to calculate median value of nodes in tree
		key_list = []
		median = None
		def h_function(node):
			if (node == None):
				return None
			key_list.append(node.key)
			h_function(node.left)
			h_function(node.right)
		h_function(self.node)
		key_list.sort()
		count = len(key_list)
		j = int(count/2)
		if count % 2:
			median =  key_list[j]
		else:
			median = (key_list[j-1] + key_list[j])/2
		return median

	def print_sum(self):
		print('Sum of elements in this subtree equals:		', self.sum_of_binary_tree())

	def print_average(self):
		print('Average of values in this subtree equals:	', self.avg_of_binary_tree())

	def print_median(self):
		print('Median of values in this subtree equals:	', self.median_of_binary_tree())

	def print_tree(self):
		#Function to recursively print out binary tree (in horizontal line - it is simplier and faster method)
		def h_function(root, count = 0):
			if (root == None):
				return
			h_function(root.right, count + 1)
			print('	'*(count) + '---| ' + str(root))
			h_function(root.left, count + 1)
		print('\nPrint out binary tree: \n')
		h_function(self.node)
		print()

#----------------------------------- FUNCTIONS -----------------------------------------------------------------
	
def go_around(root):
	#Function to go around whole tree - for testing
	def h_function(node):
		if (node == None):
			return None
		prop = Node_Prop(node)
		prop.print_tree()
		prop.print_sum()
		prop.print_average()
		prop.print_median()
		print('_' *60)

		h_function(node.left)
		h_function(node.right)
	h_function(root)
#-------------------------------------- MAIN -----------------------------------------------------------


#Lets create a tree
root = Node(5,Node(3,Node(2),Node(5)),Node(7,Node(1),Node(0,Node(2),Node(8,None,Node(5)))))
#Lets print out the tree

go_around(root)
author()
