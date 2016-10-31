from __future__ import print_function

# https://en.wikipedia.org/wiki/Red%E2%80%93black_tree
# Invariants:
# 1. A node is either red or black.
# 2. The root is black.
# 3. All leaves are black.
# 4. If a node is red, then both its children are black.
# 5. Every path from a given node to a leaf node has the same no. of black nodes.
# Results in:
# 6. Longest root->leaf path is no more than twice as long as shortest root->leaf
#    (i.e. roughly balanced)
#
# Black depth of the tree is the number of back nodes from root to any leaf
# Longest path is 2*B-1 nodes where B is the black depth of the tree
# Shortest path is B nodes

class pyRBT:
  class RBLeaf:
    def __init__(self): self.size = 0
    def isblack(self): return True
    def isred(self): return False
    def isleaf(self): return True
    def __str__(self): return "RBLeaf"
    def __len__(self): return 0
    def treestr(self): return "."

  class RBNode:
    def __init__(self,value,black=True):
      self.value = value
      self.black = black
      self.size = 1
      self.l = self.r = pyRBT.leaf
    def isblack(self): return self.black
    def isred(self): return not self.black
    def isleaf(self): return False
    def __len__(self): return self.size
    def __str__(self):
      return "RBNode("+str(self.value)+","+("black" if self.black else "red")+")"
    def treestr(self):
      col = "B" if self.black else "R"
      return "("+self.l.treestr()+","+str(self.value)+":"+col+","+self.r.treestr()+")"

  class RBIterator:
    def __init__(self,tree,reverse=False,retpaths=False):
      self.tree = tree
      self.path = []
      self.fwd = not reverse
      self.retpaths = retpaths
    def __iter__(self): return self
    def next(self): return self.__next__()
    def __next__(self):
      p = self.path
      if len(p) == 0:
        if self.tree.root.isleaf(): raise StopIteration() # empty tree
        p.append(self.tree.root)
        if self.fwd:
          while not p[-1].l.isleaf(): p.append(p[-1].l)
        else:
          while not p[-1].r.isleaf(): p.append(p[-1].r)
      elif self.fwd and not p[-1].r.isleaf():
        # Take the secondary fork (right fork when forward)
        p.append(p[-1].r)
        while not p[-1].l.isleaf(): p.append(p[-1].l)
      elif not self.fwd and not p[-1].l.isleaf():
        # Take the secondary fork (left fork when reverse)
        p.append(p[-1].l)
        while not p[-1].r.isleaf(): p.append(p[-1].r)
      else:
        # Back up the tree
        last = p.pop() # remove last leaf
        # find first parent node that used left node
        if self.fwd:
          while len(p) > 0 and last == p[-1].r: last = p.pop()
        else:
          while len(p) > 0 and last == p[-1].l: last = p.pop()
        if len(p) == 0: raise StopIteration()
      return p if self.retpaths else p[-1].value

  leaf = RBLeaf()

  def __init__(self):
    self.root = pyRBT.leaf

  def __len__(self):
    return self.root.size

  # Editing the tree voids any iterators! Do not edit the tree whilst iterating.
  def __iter__(self):
    return pyRBT.RBIterator(self,False,False)

  # Get a reverse iterator by overriding reversed(...)
  def __reversed__(self):
    return pyRBT.RBIterator(self,True,False)

  # Iterator that returns paths to each node
  def paths(self,reverse=False):
    return pyRBT.RBIterator(self,reverse,True)

  # Get a string representation of the tree
  def __str__(self):
    return self.root.treestr()

  # override the square bracket operator [] to get a value by index
  def __getitem__(self,i):
    return self.get(i)

  # setitem not defined since we don't map from key -> value
  # def __setitem__(self,key,value):

  # override the del operation to delete a value by index
  def __delitem__(self,i):
    self.pop(i)

  def __contains__(self,item):
    return self.find(item) is not None

  def clear(self):
    self.root = pyRBT.leaf

  # compare two Red Black Trees lexicographically
  # [1] < [2] < [1,1] < [1,2] < [1,2,0]
  def __cmp__(x,y):
    if len(x) != len(y): return len(x) - len(y)
    for (a,b) in zip(x,y):
      if a != b: return a-b
    return 0

  def __gt__(x,y): return x.__cmp__(y)  > 0
  def __ge__(x,y): return x.__cmp__(y) >= 0
  def __eq__(x,y): return x.__cmp__(y) == 0
  def __ne__(x,y): return x.__cmp__(y) != 0
  def __le__(x,y): return x.__cmp__(y) <= 0
  def __lt__(x,y): return x.__cmp__(y)  < 0

  # p is a path to a node path[-1]
  @staticmethod
  def _parent(path):
    return path[-2] if len(path) >= 2 else None

  @staticmethod
  def _grandparent(path):
    return path[-3] if len(path) >= 3 else None

  @staticmethod
  def _uncle(path):
    gp = pyRBT._grandparent(path)
    p = pyRBT._parent(path)
    if gp is None: return None
    return gp.r if p == gp.l else gp.l

  @staticmethod
  def _sibling(path):
    if len(path) < 2: return None
    (pa,nd) = (path[-2],path[-1])
    return pa.r if pa.l == nd else pa.l

  def _replace_node(self,pa,ch,newch):
    if pa is None: self.root = newch
    elif pa.l == ch: pa.l = newch
    else: pa.r = newch

  #
  #   pa    ->    ch
  #  /  \        /  \
  #     ch      pa
  #    /  \    /  \
  #
  # input path to parent node, upon return path points to child (new parent)
  def _rotate_left(self,path):
    pa = path.pop()
    ch = pa.r
    (pa.r, ch.l) = (ch.l, pa)
    pa.size = len(pa.l) + 1 + len(pa.r)
    ch.size = len(ch.l) + 1 + len(ch.r)
    self._replace_node(path[-1] if len(path) > 0 else None, pa, ch)
    path.append(ch)

  #
  #       pa    ->   ch
  #      /  \       /  \
  #     ch             pa
  #    / \             / \
  #
  # input path to parent node, upon return path points to child (new parent)
  def _rotate_right(self,path):
    pa = path.pop()
    ch = pa.l
    (pa.l, ch.r) = (ch.r, pa)
    pa.size = len(pa.l) + 1 + len(pa.r)
    ch.size = len(ch.l) + 1 + len(ch.r)
    self._replace_node(path[-1] if len(path) > 0 else None, pa, ch)
    path.append(ch)

  def _insert_case1(self,path):
    # print("  tree:",self)
    # print("  _insert_case1:",','.join([str(x) for x in path]))
    if len(path) == 1:
      self.root = path[0]
      self.root.black = True
    elif pyRBT._parent(path).isred():
      self._insert_case3(path)

  def _insert_case3(self,path):
    # print("  tree:",self)
    # print("  _insert_case3:",','.join([str(x) for x in path]))
    assert len(path) > 1 and path[-2].isred()
    # Assumption: parent exists and is red
    #  => therefore grandparent also exists and is black
    gp = pyRBT._grandparent(path)
    pa = pyRBT._parent(path)
    un = pyRBT._uncle(path)
    if un is not None and un.isred():
      pa.black = True
      un.black = True
      gp.black = False
      path.pop()
      path.pop()
      self._insert_case1(path) # gp is now red, deal with it
    else:
      self._insert_case4(path)

  def _insert_case4(self,path):
    # print("  tree:",self)
    # print("  _insert_case4:",','.join([str(x) for x in path]))
    gp = pyRBT._grandparent(path)
    pa = pyRBT._parent(path)
    nd = path.pop() # pop to rotate from parent
    if nd == pa.r and pa == gp.l:
      self._rotate_left(path)
      path.append(nd.l)
    elif nd == pa.l and pa == gp.r:
      self._rotate_right(path)
      path.append(nd.r)
    else:
      path.append(nd)
    self._insert_case5(path)

  def _insert_case5(self,path):
    # print("  tree:",self)
    # print("  _insert_case5:",','.join([str(x) for x in path]))
    gp = pyRBT._grandparent(path)
    pa = pyRBT._parent(path)
    nd = path[-1]
    gp.black = False
    pa.black = True
    path.pop() # pop to rotate grandparent
    path.pop()
    if nd == pa.l: self._rotate_right(path)
    else: self._rotate_left(path)

  # multiset = True allows multiple insertions of the same value
  def insert(self,item,multiset=False):
    if len(self) == 0:
      self.root = pyRBT.RBNode(item)
    else:
      p,v = [],self.root
      while not v.isleaf():
        p.append(v)
        if not multiset and item == v.value:
          v.value = item
          return
        v = (v.l if item < v.value else v.r)
      v = pyRBT.RBNode(item,black=False)
      if item < p[-1].value: p[-1].l = v
      else: p[-1].r = v
      # Need to update size
      for node in p: node.size += 1
      p.append(v)
      self._insert_case1(p)

  def extend(self,l):
    for x in l: self.insert(x)

  # remove element from a given index
  def pop(self,i=None):
    if i is None: i = len(self)-1
    p = self.getpath(i)
    return self._delete_path(p)

  # remove a given item
  def remove(self,item):
    p = self.findpath(item)
    if len(p) == 0 or p[-1].value != item:
      raise KeyError("RBT key '"+str(item)+"' not found")
    return self._delete_path(p)

  def _delete_path(self,p):
    nd = p[-1]
    val = nd.value # value that is being deleted
    # go right since we use the < and >= relations for left/right leaves
    if not nd.r.isleaf():
      v = nd.r
      p.append(v)
      while not v.l.isleaf():
        v = v.l
        p.append(v)
    elif not nd.l.isleaf():
      v = nd.l
      p.append(v)
      while not v.r.isleaf():
        v = v.r
        p.append(v)
    nd.value = p[-1].value
    for v in p: v.size -= 1
    self._delete_one_child(p)
    return val

  def _delete_one_child(self,path):
    node = path.pop()
    child = (node.l if node.r.isleaf() else node.r)
    self._replace_node(path[-1] if len(path) > 0 else None, node, child)
    path.append(child) # may be appending a leaf node, this is OK in deletion
    if node.isblack():
      if child.isred(): child.black = True
      else: self._delete_case2(path)
    # `node` is no longer in the tree

  # assume we have a parent
  def _delete_case2(self,path):
    if len(path) < 2: return
    (pa,nd,sb) = (path[-2],path[-1],pyRBT._sibling(path))
    if sb.isred():
      pa.black = False
      sb.black = True
      path.pop() # pop to rotate parent
      if nd == pa.l: self._rotate_left(path)
      else: self._rotate_right(path)
      path.append(pa)
      path.append(nd)
    self._delete_case3(path)

  def _delete_case3(self,path):
    (pa,nd,sb) = (path[-2],path[-1],pyRBT._sibling(path))
    if pa.isblack() and sb.isblack() and sb.l.isblack() and sb.r.isblack():
      sb.black = False
      path.pop()
      self._delete_case2(path) # parent
    else:
      self._delete_case4(path) # node

  def _delete_case4(self,path):
    (pa,nd,sb) = (path[-2],path[-1],pyRBT._sibling(path))
    if pa.isred() and sb.isblack() and sb.l.isblack() and sb.r.isblack():
      sb.black = False
      pa.black = True
    else:
      self._delete_case5(path)

  def _delete_case5(self,path):
    (pa,nd,sb) = (path[-2],path[-1],pyRBT._sibling(path))
    if sb.isblack():
      path.pop() # pop to rotate sibling
      path.append(sb)
      if nd == pa.l and sb.r.isblack() and sb.l.isred():
        sb.black = False
        sb.l.black = True
        self._rotate_right(path)
      elif nd == pa.r and sb.l.isblack() and sb.r.isred():
        sb.black = False
        sb.r.black = True
        self._rotate_left(path)
      path.pop() # remove sibling, re-add node
      path.append(nd)
    self._delete_case6(path)

  def _delete_case6(self,path):
    (pa,nd,sb) = (path[-2],path[-1],pyRBT._sibling(path))
    sb.black = pa.black
    pa.black = True
    path.pop() # rotate parent
    if nd == pa.l:
      sb.r.black = True
      self._rotate_left(path)
    else:
      assert nd == pa.r
      sb.l.black = True
      self._rotate_right(path)

  def find(self,item):
    v = self.root
    while not v.isleaf():
      if item == v.value: return item
      v = (v.l if item < v.value else v.r)
    return None

  # path = tree.findpath(x)
  # =>
  #   x == path[-1].value OR
  #   len(path) <= 1 OR
  #   x between path[-1].value and path[-2].value
  def findpath(self,item,p=None):
    if p is None: p,v = [],self.root
    else: v = p.pop()
    while not v.isleaf():
      p.append(v)
      if item == v.value: return p
      v = (v.l if item < v.value else v.r)
    return p

  # fetch via index
  # index is within `start` if passed
  def get(self,i,start=None):
    if start is None: v = self.root
    if i < 0: i += len(v) # allow negative indices
    if i < 0 or i >= len(v):
      raise IndexError("index out of range (%d vs 0..%d)" % (i, len(v)))
    while not v.isleaf():
      if i < len(v.l): v = v.l
      elif i == len(v.l): return v.value
      else:
        i -= len(v.l) + 1
        v = v.r
    raise RuntimeError("Internal pyRBT error")

  # fetch path via index
  # index is within `p[-1]` if passed
  def getpath(self,i,p=None):
    if p is None: p,v = [],self.root
    else: v = p[-1]
    if i < 0: i += len(v) # allow negative indices
    if i < 0 or i >= len(v):
      raise IndexError("index out of range (%d vs 0..%d)" % (i, len(v)))
    while not v.isleaf():
      p.append(v)
      if i < len(v.l): v = v.l
      elif i == len(v.l): return p
      else:
        i -= len(v.l) + 1
        v = v.r
    raise RuntimeError("Internal pyRBT error")

  # Get the index of an item
  def index(self,item,start=None):
    if start is None: v = self.root
    i = 0
    idx = None
    while not v.isleaf():
      if item < v.value: v = v.l
      elif item == v.value:
        # found one instance, look for earlier ones
        idx = i+len(v.l)
        v = v.l
      else:
        i += len(v.l) + 1
        v = v.r
    if idx is None: raise KeyError('Key not found: '+str(item))
    return idx

  # Check data structure integrity by checking invariants are met
  def check(self):
    assert (len(self) == 0) == self.root.isleaf() # size is zero only if empty
    assert self.root.isblack() # root node is black
    nblack = -1
    npaths = 0
    for p in self.paths():
      # print("Check:",'->'.join([str(x) for x in p]))
      assert not p[-1].isleaf() or p[-1].isblack() # all leaf nodes are black
      if p[-1].isred():
        # all red nodes have only black children
        assert p[-1].l.isblack() and p[-1].r.isblack()
      # Every path from the the root has the same number of black nodes
      if p[-1].l.isleaf() or p[-1].r.isleaf():
        ntmpb = sum([ x.isblack() for x in p ]) + 1
        assert nblack == -1 or nblack == ntmpb
        nblack = ntmpb
      npaths += 1
    assert(npaths == len(self))
    # print('nblack:',nblack,'npaths:',npaths)


import random

def _test_rbt(nums):
  tree = pyRBT()
  vals = [] # sorted values in the tree
  # insert values into the tree
  for v in nums:
    tree.insert(v,True)
    vals.append(v)
    vals.sort()
    tree.check()
    assert sum([ x==y for x,y in zip(vals,iter(tree)) ]) == len(vals)
    # Test indexing
    for (i,v) in enumerate(vals): assert(tree[i] == v)
    for (i,v) in enumerate(vals): assert tree.index(v) == vals.index(v)
  # remove the values in a random order
  rvals = list(nums)
  random.shuffle(rvals)
  for v in rvals:
    tree.remove(v)
    vals.remove(v)
    tree.check()
    assert sum([ x==y for x,y in zip(vals,iter(tree)) ]) == len(vals)
    # Test indexing
    for (i,v) in enumerate(vals): assert tree[i] == v
    for (i,v) in enumerate(vals): assert tree.index(v) == vals.index(v)
  # Re-build sorted remove forwards
  sortedvals = sorted(vals)
  tree.extend(sortedvals)
  for v in sortedvals: assert sortedvals.remove(v) == tree.remove(v)
  assert len(tree) == 0
  tree.check()
  # Re-build sorted remove backwards
  sortedvals = sorted(vals)
  tree.extend(sortedvals)
  for v in reversed(sortedvals): assert sortedvals.remove(v) == tree.remove(v)
  assert len(tree) == 0
  tree.check()

def _test_rbt_comparison():
  print("Testing RBT comparison...")
  abc = pyRBT()
  xyz = pyRBT()
  abc.extend([1,2,3,9])
  xyz.extend([9,3,2,1])
  assert abc == xyz and abc >= xyz and abc <= xyz
  assert not (abc > xyz) and not (abc < xyz) and not (abc != xyz)
  abc.remove(3)
  # abc < xyz
  assert abc < xyz and xyz > abc and abc != xyz
  assert not (abc > xyz) and not (abc >= xyz) and not (abc == xyz)
  assert not (xyz < abc) and not (xyz <= abc)
  abc.clear()
  xyz.clear()
  # abc == xyz
  assert abc == xyz and abc <= xyz and abc >= xyz
  assert not (abc != xyz) and not (abc < xyz) and not (abc > xyz)
  # abc > xyz
  abc.insert(1)
  assert abc > xyz and xyz < abc and abc != xyz
  assert not (abc < xyz) and not (abc <= xyz) and not (abc == xyz)
  assert not (xyz > abc) and not (xyz >= abc)


def _test_red_black_tree():
  print("Testing RBT")
  tree = pyRBT()

  # Insert [1,2,...,N]
  print("Doing mixed automated tests...")
  vals = list(range(1,100))
  _test_rbt(vals) # 1..N in sorted order
  random.shuffle(vals)
  _test_rbt(vals) # 1..N in random order
  _test_rbt([]) # test epmty set
  _test_rbt([random.randrange(100) for x in range(200)]) # random numbers with repeats
  _test_rbt([random.randrange(1000) for x in range(200)]) # random numbers with repeats
  _test_rbt([1]*10) # multiple 1s

  # Test python features
  vals = list(range(1,10))
  sortedvals = sorted(list(vals)) # sorted copy
  random.shuffle(vals)
  for v in vals:
    tree.insert(v,True)
    tree.check()
  # -- Testing iterate paths --
  # print("Testing path iteration...")
  # for path in tree.paths():
  #   print(' ','->'.join([str(x) for x in path]))
  # print("Testing path iterating reversed...")
  # for path in tree.paths(reverse=True):
  #   print(' ','->'.join([str(x) for x in path]))
  print("Testing value iteration...")
  assert ','.join([str(x) for x in tree]) == ','.join([str(x) for x in sortedvals])
  print("Testing value iteration reversed...")
  assert ','.join([str(x) for x in reversed(tree)]) == ','.join([str(x) for x in reversed(sortedvals)])
  print("Testing tree[i]...")
  for i in range(len(tree)): assert tree[i] == sortedvals[i]
  print("Testing find...")
  for f in [-1,0.5,len(vals)+1,len(vals)-3.6]: assert tree.find(f) is None
  print("Testing tree.index(i)...")
  for (i,v) in enumerate(sortedvals): assert tree.index(v) == sortedvals.index(v)
  print("Checking tree...")
  tree.check()
  # Test removing random nodes with pop
  while len(sortedvals) > 0:
    i = random.randrange(len(sortedvals))
    assert tree.pop(i) == sortedvals.pop(i)
    tree.check()
  # re-build tree and sorted list
  tree.extend(vals)
  sortedvals = sorted(list(vals)) # sorted copy
  assert sum([ x == y for (x,y) in zip(iter(tree),sortedvals)]) == len(vals)
  print("Testing remove...")
  for v in vals:
    tree.remove(v)
    tree.check()
  _test_rbt_comparison()
  print("Looks like the tests all passed...")

if __name__ == '__main__':
  _test_red_black_tree()

del _test_rbt
del _test_red_black_tree
