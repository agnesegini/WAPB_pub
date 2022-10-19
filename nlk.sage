##computing weightwise nonlinearity on the slices

#compute the hamming weight of a vector x
def hw(x):
     x=vector(ZZ,x)
     w=vector(ones_matrix(ZZ,1,len(x)))
     return w*x
     
#returns the hamming distance of the F2's vectors c and r  
def h_dist(c,r):
   return hw(c+r)


#given as input two integers n and k retunrs the linear code C over F2, that is a spherically punctured Reed Muller code of order 1 of lenght v=(n choose k)
def Pkn(k,n):
  b=binomial(n,k)
  M=ones_matrix(GF(2),n+1,b) 
  E=Ekn(k,n)
  for i in range(b):
    M[1:,i]=int_bin_n(E[i],n) ##see supp_prop.sage
  C=LinearCode(M)
  return C
  
#return the list of values of the distance between the vector r and any vector of the code C
def dist_all(C,r):
  M=C.generator_matrix()
  d=C.dimension()
  it=prodit(range(2), repeat=d)
  return map(lambda x: h_dist(vector(GF(2),x)*M,r) ,it)  
 
#returns the restriction to S of the Boolean function f 
def restriction(f,S):
    fT=f.truth_table('int')
    return vector(GF(2),[fT[i] for i in S])


#return the distance between the binary vectors x*M and t
def fun_par(M,r,x):
  return h_dist(vector(GF(2),x)*M,r)      


#return the list of values of the distance between the vector r and any vector of the code C
#this runs in parallel         
def dist_all_parallel(C,r):
  pool = Pool(cpu)
  M1=C.generator_matrix()
  d=C.dimension()
  it=prodit(range(2), repeat=d)
  fun=partial(fun_par,M1,r)
  D=pool.imap(fun,it)#,   map(lambda x: vector(GF(2),x)*M ,it) 
  pool.close()
  pool.join()
  return list(D)

#returns the weightwise nonlinearity of n-variable f on the slice Ekn
#via computing the distance bewteen the spherically punctured Reed Muller code of order 1 of lenght v=(n choose k)  and the support of the restriction of f over the slide Ekn  as suggested in  https://ia.cr/2022/408
#this is sequential

def NLk(k,f):
    n=f.nvariables()
    if k<0 or k>n: return 0
    vf=restriction(f,Ekn(k,n))
    return min(dist_all(Pkn(k,n),vf))
  
#returns the weightwise nonlinearity of n-variable f on the slice Ekn  
#via computing the distance bewteen the spherically punctured Reed Muller code of order 1 of lenght v=(n choose k)  and the support of the restriction of f over the slide Ekn as suggested in https://ia.cr/2022/408
#this runs in parallel

def NLk_par(k,f):
    n=f.nvariables()
    if k<0 or k>n: return 0
    vf=restriction(f,Ekn(k,n))
    return min(dist_all_parallel(Pkn(k,n),vf))  

#returns the weightwise nonlinearity of n-variable f on the slice Ekn  
#via walsh tranform
#his is sequential  
def NLk_w(k,f):
    n=f.nvariables()
    b=binomial(n,k)
    w=max([abs(walsh(f,a,Ekn(k,n))) for a in range(2^n)])
    return b/2-w/2
  
