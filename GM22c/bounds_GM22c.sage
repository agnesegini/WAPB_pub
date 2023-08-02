"""
  Computation bounds on the nonlinearity of a WPB function
"""
def _K(k,x,n):
  return codes.bounds.krawtchouk(n,2, k,x)

def K(k,x,n):
  """Return the Krawtchouk Polynomial K_k(x, n) (Definition 10 [GM22c])

  Args:
      k (int): degree
      x (int): evaluation point
      n (int): 

  Returns:
      _type_: value of K_k(x, n) 
  """  
  return sum([ (-1)^j*binomial(x,j)*binomial(n-x,k-j) for j in range(k+1)])
  
def mu_kn(k,m):
  """Return the minimum weightwise nonlinearity of WPB function on a specific slice mu_kn (Definition 9 [GM22a])

  Args:
      k (int): minimum weightwise nonlinearity on the slice k
      m (int): 2^m is the number of variables

  Returns:
      int: mu_kn
  """  
  n=2^m
  L=[abs(K(k,l+1,n))/2 for l in range(n/2)]
  return min(L)
   
def min_list_K(t,n):
  return min([abs(K(2*t,l+1,n)) for l in range(n/2)])
  
def Tl(l,n):
  T=abs(K(n//2,l,n))/2
  for k in range(1,n//2): T+=abs(K(k,l,n))
  return T
  
def Bm(m):
  """Lower bound on the nonlinearity of WPB functions (Theorem 1 [GM22c])

  Args:
      m (int): 2^m is the number of variables

  Returns:
      int: lower bound B_m
  """   
  n=2^m
  T=Tl(0,n)
  mi=T+1
  for l in range(1,n//2+1):
    T= Tl(l,n)
    mi=min([mi,T+1+(l%2),T+((l+1)%2)]) 
  return mi
  
