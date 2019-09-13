
tname='_yeast'

tree=$(<yeast_tree.new)

tchr=$(<yeast_nodes.txt)
readarray atchr < yeast_nodes_split.txt

lchr=$(<yeast_leafs.txt)

step=0.001


