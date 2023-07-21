#!/usr/bin/python3
## SageMath version 9.5

"""
SLICES,WPB,WAPB,SWAPB
"""

def Ekn(k,n):
  """Retuns integers between 0 and 2^n whose binary rapresentation has hamming weight k.

  Args:
      k (int): Hamming weight
      n (int): bitsize

  Returns:
      list : integers between 0 and 2^n with Hw k.
  """  
  return list(filter(lambda d : (hw(bin(d)[2:])==k), range(2^n)))

 # 
def Ikn(k,n):
  """Returns the slice indicator function (phi_kn) 

  Args:
      k (int): Hamming weight
      n (int): bitsize

  Returns:
      BooleanFunction : slice indicator function
  """
  L=2^n*[0]
  for i in Ekn(k,n): L[i]=1
  return BooleanFunction(L)

def a_Ikn(k,n):
  """Returns the algebrain normal form of the slice indicator functions phi_kn"""
  L=2^n*[0]
  for i in Ekn(k,n): L[i]=1
  return BooleanFunction(L).algebraic_normal_form()
   
def least_bit(x):
  """Returns  the position least non-zero bit of an integer """
  c=x.binary()
  l=len(c)
  i=0
  while True:
    if '0'==c[l-1-i]: i+=1 
    else: return i 


def supp(f):
  """Returns the support of the Boolean function f, ie the list of integers a such that f(a)=1

  Args:
      f (BooleanFunction): 

  Returns:
      list : support of f
  """  
  n=f.nvariables()
  fT=f.truth_table(format='int')
  return list(filter(lambda x : (fT[x]==1), range(2^n)))
  #[x for x in xrange(2^n) if fT[x]==1 ]
  
def size_supp(f):
  """Returns cardinality of the support of the Boolean function f  """
  return hw(f.truth_table(format='int'))

# 
def suppk(f,k,iter=False):
  """Computes the support of the Boolean function f restricted to the slice Ekn 
  Args:
      f (Boolean function): _description_
      k (_type_): _description_
      iter (bool, optional): returns an iterator if True. Defaults to False.

  Returns:
      list or iter : support of the Boolean function f restricted to the slice Ekn
  """ 
  n=f.nvariables()
  #E=Ekn(k,n)
  fT=f.truth_table(format='int')
  if iter: return filter(lambda x : (fT[x]==1), Ekn(k,n))
  return list(filter(lambda x : (fT[x]==1), Ekn(k,n)))
 
def size_suppk(f,k):
  """Returns cardinality of the support of the Boolean function f restricted to the slice Ekn """
  return len(suppk(f,k))  

#### Verification Properties 

def int_bin_n(d,n):
  """Given an integer d=sum d_{n-i} 2^i return the vector (d_0,...,d_{n-1})"""
  d=bin(d)[2:]
  ld=len(d)
  assert (n-ld)>=0, ld
  d=(n-ld)*'0'+d
  return vector(ZZ,d)

def __Lin(a,n):
  va=vector(GF(2),int_bin_n(a,n))
  return BooleanFunction([va.dot_product(vector(GF(2),int_bin_n(d,n))) for d in range(2^n)])


def walsh_tr_0(f,S=0):
  """Returns Walsh transform  of f restricted to S at zero, if S=0 this is the Walsh transform at zero"""
  n=f.nvariables()
  if not S: S=range(2^n)
  w=0
  fT=f.truth_table(format='int')
  for i in S:
    w+=(-1)^(fT[i])
  return w


def walsh(f,a,S=0):
  """Returns Walsh transform of f restricted to S at a, if S=0 this is the Walsh transform at a"""
  n=f.nvariables()
  fa=f+__Lin(a,n)
  if not S: S=range(2^n)
  w=0
  fT=fa.truth_table(format='int')
  for i in S:
    w+=(-1)^(fT[i])
  return w
  

def is_SWAPB(f):
  """ Checks if a function is SWAPB. The verification done via  Walsh transform.
  Args:
      f (BooleanFunction)

  Returns:
      bool : True iff f is SWAPB
  """  
  n=f.nvariables()
  for k in range(0,n//2):
    b=binomial(n,k)
    if b%2==0 and walsh_tr_0(f,Ekn(k,n))!=0: False,k
    if b%2 and walsh_tr_0(f,Ekn(k,n))!=1: False,k
  for k in range(n//2,n+1):
    b=binomial(n,k)
    if b%2==0 and walsh_tr_0(f,Ekn(k,n))!=0: False,k
    if b%2 and walsh_tr_0(f,Ekn(k,n))!=-1: False,k
  return True

def is_WAPB(f):
  """ Checks if a function is WAPB. The verification done via  Walsh transform.
  Args:
      f (BooleanFunction)

  Returns:
      bool : True iff f is WAPB
  """  
  n=f.nvariables()
  for k in range(1,n/2):
    b=binomial(n,k)
    w=walsh_tr_0(f,Ekn(k,n))
    if b%2==0 and w!=0: False,k
    if b%2 and not (w in [-1,1]): False,k
  return True  
      
 
def is_WPB(f,verbose=False):
  """ Checks if a function is WAPB. The verification done via supports computation.
  Args:
      f (BooleanFunction)

  Returns:
      bool: True iff f is WPB
      (bool,int): (False, k) k is a slice cause failure
  """  
  n=f.nvariables()
  assert n%2==0 and n.is_prime_power(), "number of variables not a power of 2!"
  fT=f.truth_table(format='int')
  if fT[0] or not fT[-1]: return False,(0,n)
  for k in range(1,n):
    if verbose: print(k)
    b=binomial(n,k)
    if size_suppk(f,k)!=b/2: return False,k
  return True,True

def is_WAPB_2(f):
  """ Checks if a function is WAPB. The verification done via supports computation.
  Args:
      f (BooleanFunction)

  Returns:
      bool: True iff f is WPB
      (bool,int): (False, k) k is a slice cause failure
  """  
  fT=f.truth_table(format='int')
  n=f.nvariables()
  for k in range(1,n):
    b=binomial(n,k)
    ssk=size_suppk(f,k)
    if not b%2 and  ssk!=b/2: return False,k
    if  b%2 and (ssk not in [(b-1)/2,(b+1)/2] ): return False,k
  return True


def is_SWAPB_2(f):
  """ Checks if a function is SWAPB. The verification done via support computation.
  Args:
      f (BooleanFunction)

  Returns:
      bool : True iff f is SWAPB
  """    
  fT=f.truth_table(format='int')
  n=f.nvariables()
  if fT[0] or not fT[-1]: return False, 0
  for k in range(1,n//2+1):
    b=binomial(n,k)
    ssk=size_suppk(f,k)
    if not b%2 and  ssk!=b/2: return False,k
    if  b%2 and (ssk not in [(b-1)/2] ): return False,k
  for k in range(n//2+1,n):
    b=binomial(n,k)
    ssk=size_suppk(f,k)
    if not b%2 and  ssk!=b/2: return False,k
    if  b%2 and (ssk not in [(b+1)/2] ): return False,k
  return True  
 


  

