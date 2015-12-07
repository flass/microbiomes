#!/usr/bin/python
# -*- coding: utf8 -*-

from HTMLParser import HTMLParser
import os
import copy

class TaxonomyNode(object):
	def __init__(self, dattr={}, name=None):
		self.__children = []
		self.__father = None
		self.__dict__.update(dattr)
		if name: self.__name = name
		
	def __getattr__(self, attr):
		if attr=='name':
			return self.__name
		#~ elif attr in self.__dict__:
			#~ return self.__dict__[attr]
		#~ else:
			#~ #raise AttributeError(attr)
			#~ return None
		else:
			return self.__dict__.get(attr)
		
	def __iter__(self):
		return self.generator()		
		
	def generator(self):
		children = self.get_all_children()
		for child in children:
			yield child
			
	def get_all_children(self):
		"""Return the list of all nodes below the Node, including itself, in a pre-order traversal"""
		a=[self]
		#~ if self.__children!=[]:
		for i in self.__children:
			a+=i.get_all_children()
		return a
		
	def go_father(self):
		return self.__father
		
	def get_children(self):
		return self.__children
		
	def _addchild(self, node):
		self.__children.append(node)
		
	def _addfather(self, node):
		self.__father = node
		
	def newchild(self, dattr={}, name=None):
		node = TaxonomyNode(dattr=dattr, name=name)
		self._addchild(node)
		node._addfather(self)
		return node
		
	def updateattr(self, dattr):		
		self.__dict__.update(dattr)
		
	def __getitem__(self, key):
		if self.name==key: return self
		for child in self.__children:
			if child.name==key:
				return child
			else:
				grandchild=child[key]	  # <==> grandchild = child.__getitem__(key)   (recursive function)
				if grandchild:
					return grandchild
		return None
		
	def get(self, key, default=None):
		g = self[key]
		if g: return g
		else: return default
		
	def taxidgetnode(self, key):
		if self.taxon==key: return self
		for child in self.__children:
			if child.taxon==key:
				return child
			else:
				grandchild=child.taxidgetnode(key)	  #    (recursive function)
				if grandchild:
					return grandchild
		return None
		
	def __str__(self):
		if self.__father: fatname = self.__father.name
		else: fatname = None
		return "TaxonomyNode '%s'\nfather: %s\n%d children: %s\nattributes: %s"%(self.name, str(fatname), len(self.__children), str([child.name for child in self.__children]), str(self.__dict__))
		
	def get_cumulative_values(self, attr, exclude=[], filterNull=True, storeAsAttr=False, applyFun=None):
		lval = []
		if storeAsAttr:
			precomp = self.__getattr__(storeAsAttr)
			if precomp: return precomp
		if not self in exclude:
			lval += [self.__getattr__(attr)]
			for child in self.__children:
				lval += child.get_cumulative_values(attr, exclude=exclude, filterNull=filterNull, storeAsAttr=storeAsAttr, applyFun=applyFun)	# recursive call
		if filterNull: lval = [val for val in lval if val]
		if applyFun: lval = [applyFun(lval)]
		if storeAsAttr: self.__setattr__(storeAsAttr, lval)
		return lval
		
	def map_to_node(self, lnodes, idattrtype='taxon'):
		"""finds deepest node/clade in tree that includes all leaves in input label list = MRCA = most recent common ancestor of the leaf set
		
		uses a desending algorithm
		"""
		si = set(lnodes)
		ll = set(self.get_cumulative_values(idattrtype))
		if not si <= ll:
			if force:
				si = si & ll
			else:
				raise IndexError, "Input leaf set is larger than leaves present in reference tree:\nsi: %r\nll: %r"%(si, ll)
		if not si:
			return None
		clade = self
		scl = ll
		#~ print "si", si
		#~ print "scl", scl
		while si <= scl:
			for child in clade.children:
				scl = child.get_cumulative_values(idattrtype)
				if si <= scl:
					#~ print "scl", scl
					clade = child
					break
			else:
				break
		return clade
		
	def get_cumulative_members(self, exclude=[]):
		#~ lmemb = []
		#~ if not self in exclude:
			#~ if self.members: lmemb += [self.members]
			#~ for child in self.__children:
				#~ lmemb += child.get_cumulative_members(exclude=exclude)	# recursive call
		lmemb = self.get_cumulative_values('members', exclude=exclude)
		return lmemb
		
	def write_attr_to_table(self, lattr, nfout, sep='\t'):
		if isinstance(nfout, file): fout = nfout
		elif os.path.exists(nfout): fout = open(nfout, 'w')
		else: raise ValueError, 'wrong file descriptor %s'%(str(nfout))
		fout.write(sep.join([str(self.get(attr) for attr in lattr)])+'\n')
		for child in self.__children: child.write_attr_to_table(lattr, fout, sep=sep)	# recursive call
		if not self.go_father(): fout.close()		

class KronaHTMLParser(HTMLParser):
	
	def __init__(self, verbose=False):
		HTMLParser.__init__(self)
	#~ def __init__(self, msg, position=(None, None)):
		#~ super(KronaHTMLParser, self).__init__(msg, position)
		self.__currtag = None
		self.__prevtag = None
		### enriched attributes towards taxonomy parsing
		self.__topnode = None	# top node of the parsed taxonomy, contains all others
		self.__currnode = None	# current taxonomic node being filled
		self.__taxonomy = []	# list of taxonomy trees, to be returned
		self.__verbose = verbose
	
	def __getattr__(self, attr):
		if attr=='topnode':
			return self.__topnode
		elif attr=='currnode':
			return self.__currnode
		elif attr=='currtag':
			return self.__currtag
		elif attr=='prevtag':
			return self.__prevtag
		elif attr=='taxonomy':
			return self.__taxonomy
		elif attr in self.__dict__:
			return self.__dict__[attr]
		else:
			raise AttributeError(attr)
	
	def __setattr__(self, attr, value):
		if attr=='topnode':
			self.__topnode = value
		elif attr=='currnode':
			self.__currnode = value
		elif attr=='currtag':
			self.__currtag = value
		elif attr=='prevtag':
			self.__prevtag = value
		else:
			self.__dict__[attr] = value
	
	def handle_starttag(self, tag, attrs):
		# duplicate attribute surveying the current tag ; 'self.lasttag' is internal and should not be touched
		self.currtag = copy.copy(tag)
		#### taxonomy parsing handler
		if self.__verbose: print 'tag', tag, 'currtag', self.currtag, 'prevtag', self.prevtag
		if tag=='node':
			dattr = dict(attrs)
			if self.__verbose: print tag, dattr
			if not self.currnode:
				# initiates the taxonomy with a root node
				taxnode = TaxonomyNode(dattr=dattr)
				if not self.topnode: self.topnode = taxnode
				else: raise IndexError, 'cannot replace top node'
			else:
				taxnode = self.currnode.newchild(dattr=dattr)
			self.currnode = taxnode
			if self.__verbose: print "# update currnode to child:"
			if self.__verbose: print self.currnode
		elif tag=='val':
			# 'val' is a bogus tag, reset it to the previous one
			if self.__verbose: print '  change currtag', self.currtag, 'to', self.prevtag
			self.currtag = copy.copy(self.prevtag)
		# end of tag parsing, update previous tag
		if self.__verbose: print 'update prevtag', self.prevtag, 'to', self.currtag
		self.prevtag = copy.copy(self.currtag)
			
	def handle_data(self, data):
		if self.currnode:
			dat = data.strip(' \n')
			if dat: self.currnode.updateattr({self.currtag:dat})
	
	def handle_endtag(self, tag):
		#### taxonomy parsing handler
		if tag=='node':
			# closes the node and goes up to its father
			self.currnode = self.currnode.go_father()
			if self.__verbose: print "# update currnode to fat:"
			if self.__verbose: print self.currnode
			if self.currnode==None:
				self.__taxonomy.append(self.topnode)
				if self.__verbose: print '# reached top of taxonomy at node:'
				if self.__verbose: print self.topnode
				if self.__verbose: print '# store it and resume'
				self.topnode = None
	
