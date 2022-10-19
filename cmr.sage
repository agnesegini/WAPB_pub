## building  CMR WAPB functions from https://ia.cr/2017/097


# recursive call for CMR contructions 
def CMR_call(s,V,R):
  if s==2: return var(V[0])
  if s%2: return  CMR_call(s-1,V,R)
  if s.is_prime_power(): 
    d=log(s,2)
    L=[R(V[s-2-i]) for i in range(2^(d-1)) ]
    return CMR_call(s-1,V,R)+R(V[s-3])+prod(L)
  d=least_bit(s)
  L=[R(V[s-2-i]) for i in range(2^(d)) ]
  return CMR_call(s-1,V,R)+var(V[s-3])+prod(L)  

#return the CMR WAPB construction in n variables
def CMR(n):
  R=BooleanPolynomialRing(names=['x'+str(i) for i in range(n)])
  V=R.variable_names()
  f=CMR_call(ZZ(n),V,R)
  return BooleanFunction(R(f))
  
