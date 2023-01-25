#!/usr/bin/python3
## SageMath version 9.5


####SLICES,WPB,WAPB,SWAPB

##supports, indicator functions

# retuns the list of intregers between 0 and 2^n whose binary rapresentation has hamming weight k. 
# This represent the list of indices of the element of the slice Ekn in range(0,2^n)
def Ekn(k,n):
  return list(filter(lambda d : (hw(bin(d)[2:])==k), range(2^n)))
  #[d for d in xrange(2^n) if hw(bin(d)[2:])==k]

 # returns the slice indicator functions phi_kn 
def Ikn(k,n):
  L=2^n*[0]
  for i in Ekn(k,n): L[i]=1
  return BooleanFunction(L)

 # returns the algebrain normal form of the slice indicator functions phi_kn 
def a_Ikn(k,n):
  L=2^n*[0]
  for i in Ekn(k,n): L[i]=1
  return BooleanFunction(L).algebraic_normal_form()
  
# returns  the position least non-zero bit of an integer 
def least_bit(x):
  c=x.binary()
  l=len(c)
  i=0
  while True:
    if '0'==c[l-1-i]: i+=1 
    else: return i 


# returns  the support of the Boolean function f, ie the list of integers a such that f(a)=1
def supp(f):
  n=f.nvariables()
  fT=f.truth_table(format='int')
  return list(filter(lambda x : (fT[x]==1), range(2^n)))
  #[x for x in xrange(2^n) if fT[x]==1 ]
  
# returns cardinality of the support of the Boolean function f 
def size_supp(f):
  return hw(f.truth_table(format='int'))

# returns  the support of the Boolean function f restricted to the slice Ekn 
def suppk(f,k):
  n=f.nvariables()
  #E=Ekn(k,n)
  fT=f.truth_table(format='int')
  return list(filter(lambda x : (fT[x]==1), Ekn(k,n)))

# returns cardinality of the support of the Boolean function f restricted to the slice Ekn  
def size_suppk(f,k):
  return len(suppk(f,k))  

##verification properties

# given an integer d=sum d_{n-i} 2^i return the vector (d_0,...,d_{n-1})
def int_bin_n(d,n):
  d=bin(d)[2:]
  ld=len(d)
  assert (n-ld)>=0, ld
  d=(n-ld)*'0'+d
  return vector(ZZ,d)

# 
def __Lin(a,n):
  va=vector(GF(2),int_bin_n(a,n))
  return BooleanFunction([va.dot_product(vector(GF(2),int_bin_n(d,n))) for d in range(2^n)])

#returns Walsh transform  of f restricted to S at zero, if S=0 this is the Walsh transform at zero
def walsh_tr_0(f,S=0):
  n=f.nvariables()
  if not S: S=range(2^n)
  w=0
  fT=f.truth_table(format='int')
  for i in S:
    w+=(-1)^(fT[i])
  #print w
  return w

#returns Walsh transform of f restricted to S at a, if S=0 this is the Walsh transform at a
def walsh(f,a,S=0):
  n=f.nvariables()
  fa=f+__Lin(a,n)
  if not S: S=range(2^n)
  w=0
  fT=fa.truth_table(format='int')
  for i in S:
    w+=(-1)^(fT[i])
  #print w
  return w
  
#returns True if f is SWAPB; False otherwise. 
#verification done via  Walsh transform
def is_SWAPB(f):
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

#returns True if f is WAPB; False otherwise. 
#verification done via  Walsh transform  
def is_WAPB(f):
  n=f.nvariables()
  for k in range(1,n/2):
    b=binomial(n,k)
    w=walsh_tr_0(f,Ekn(k,n))
    #print w
    if b%2==0 and w!=0: False,k
    if b%2 and not (w in [-1,1]): False,k
  return True  
      

#returns True if f is WPB; False otherwise. 
#verification done via supports computation
def is_WPB(f,verbose=False):
  n=f.nvariables()
  assert n%2==0 and n.is_prime_power(), "not a power of 2!"
  fT=f.truth_table(format='int')
  if fT[0] or not fT[-1]: return False,(0,n)
  for k in range(1,n):
    if verbose: print(k)
    b=binomial(n,k)
    if size_suppk(f,k)!=b/2: return False,k
  return True,True

#returns True if f is WAPB; False otherwise. 
#verification done via supports computation
def is_WAPB_2(f):
  fT=f.truth_table(format='int')
  n=f.nvariables()
  for k in range(1,n):
    b=binomial(n,k)
    ssk=size_suppk(f,k)
    if not b%2 and  ssk!=b/2: return False,k
    if  b%2 and (ssk not in [(b-1)/2,(b+1)/2] ): return False,k
  return True

#returns True if f is SWA.PB; False otherwise. 
#verification done via supports computation
def is_SWAPB_2(f):
  fT=f.truth_table(format='int')
  n=f.nvariables()
  if fT[0] or not fT[-1]: return False
  for k in range(1,n//2+1):
    b=binomial(n,k)
    ssk=size_suppk(f,k)
    if not b%2 and  ssk!=b/2: return False
    if  b%2 and (ssk not in [(b-1)/2] ): return False,k
  for k in range(n//2+1,n):
    b=binomial(n,k)
    ssk=size_suppk(f,k)
    if not b%2 and  ssk!=b/2: return False,k
    if  b%2 and (ssk not in [(b+1)/2] ): return False,k
  return True  
 


  

