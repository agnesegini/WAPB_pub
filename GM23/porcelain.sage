

def kaolinite(n): 
   R=BooleanPolynomialRing(names=['x'+str(i) for i in range(n)])
   f=R('x0*(x1+x2)')
   return BooleanFunction(f)
   
def porcelain_0(n):
   F=kaolinite(n)
   G=~F
   N=list(G.truth_table('int'))
   N[0]=0
   N[-1]=1
   wpb, t =is_WPB(BooleanFunction(N))
   #L=[]
   while not wpb :
      print(t)
      E=Ekn(t,n)
      o=walsh_tr_0(G,E)//2
      print(t,o)
      S=suppk(G,t)
      #L+=[(t,o,len(CS), binomial(len(CS),o))]
      if o>0: 
         CS = [x for x in E if x not in S]
         for i in range(o):
            N[CS[i]]=1
      else: 
         for i in range(-o):
            N[S[i]]=0
      y=BooleanFunction(N)
      wpb, t =is_WPB(y) 
   nly=y.nonlinearity()
   return y,nly,wpb
  #print(L)
  
def porcelain_rand(n, flipped=False):
   F=kaolinite(n)
   p=construction1_GM22c(~F)
   if flipped: return p
   else: return p[0]  
   
   
def walsh_all_slices(f):
  n=f.nvariables()
  return [walsh_tr_0(f,Ekn(i,n)) for i in range(n+1)]

def numPorcelain(n):
   num=binomial(n,n/2)
   b=binomial(n,2)
   num*=binomial(b-2,b/2-2)
   for k in range(3,n):
      b=binomial(n,k)
      ob=binomial(n-3,k-2)
      num*=binomial(b-2*ob,b//2-2*ob)
   return num
   

