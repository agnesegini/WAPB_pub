#load("main.sage")


def _K(k,x,n):
  return codes.bounds.krawtchouk(n,2, k,x)

def K(k,x,n):
  return sum([ (-1)^j*binomial(x,j)*binomial(n-x,k-j) for j in range(k+1)])
  
def mu_kn(k,m):
  n=2^m
  L=[abs(K(k,l+1,n))/2 for l in xrange(n/2)]
  #print L
  return min(L)
   
def min_list_K(t,n):
  return min([abs(K(2*t,l+1,n)) for l in xrange(n/2)])
  
def Tl(l,n):
  T=abs(K(n//2,l,n))/2
  for k in range(1,n//2): T+=abs(K(k,l,n))
  return T
  
def Bm(m): 
  n=2^m
  T=Tl(0,n)
  mi=T+1
  for l in range(1,n//2+1):
    T= Tl(l,n)
    mi=min([mi,T+1+(l%2),T+((l+1)%2)]) 
    #print(mi)
  return mi
  
