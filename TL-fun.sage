from random import sample
#c=[1,0,0,0,-1]
def LT_Wi(k,i,rel):
   n=2*k
   Fk1=itertools.product(range(2),repeat=k)
   if rel=="Eq":
      W=[]
      for x in Fk1:
         Fk2=itertools.product(range(2),repeat=k)
         hx=hw(x)
         if hx==i//2:
           for y in Fk2:
             hy=hw(y)
             if (hy==i//2) and (hx+hy==i): W.append(x+y) 
   elif rel=="Gr":
      W=[]
      for x in Fk1:
         Fk2=itertools.product(range(2),repeat=k)
         hx=hw(x)
         for y in Fk2:
             hy=hw(y)
             if (hx+hy==i) and hx>hy:
               W.append(x+y)
   return W
 
def LT_WG(k):
   n=2*k 
   W=set(LT_Wi(k,0,"Gr"))
   for i in range(0,n+1):
      Wi=LT_Wi(k,i,"Gr")
      for w in Wi: W.add(w)
   return W
       
def LT_Ui(k,i,c):
   if mod(i,2)==1: return [] 
   else: 
      W=LT_W(k,i,"Eq")
      w=len(W)
      ui=(w+c)//2
      U=random.sample(W, ui)
      assert len(U)==ui
      return U

def LT_U(k,c):
   n=2*k 
   U=set(LT_Ui(k,0,c[0]))
   for i in range(1,k+1):
      U2i=LT_Ui(k,2*i,c[i])
      for w in U2i: U.add(w)
   return U

def LT_is_valid_seq(c,k):
   sc=0
   for j in range(k+1):
      bi=binomial(k,j)
      #print(bi,c[j])
      if abs(c[j])!=mod(bi,2):
         return False
      sc+=c[j] 
   if sc!=0: return False
   return True

def LT_from_seq(c,k):
   n=2*k
   assert LT_is_valid_seq(c,k)
   WG=LT_WG(k)
   U=LT_U(k,c)
   U.union(WG)
   print(len(U))
   Lf=[]
   for d in range(2^n):
      v=bin(d)[2:]
      v=(n-len(v))*"0"+v
      v=tuple(vector(v,ZZ))
      if v in U: Lf+=[1]
      else: Lf+=[0]
   print(Lf)
   return BooleanFunction(Lf)
   
def LT_gr(v,k):
   x=v[:k]
   y=v[k:]
   return 
   hw(x)>hw(y)

def LT_bin_to_int(v): return int("".join(str(i) for i in v),2)
   
def LT_WPB_rand(n):
   k=n//2
   vf=vector(ZZ,2**n)
   vf[-1]=1
   for d in range(2^n):
      v=bin(d)[2:]
      v=(n-len(v))*"0"+v
      x=v[:k]
      y=v[k:]
      hx=hw(x)
      hy=hw(y)
      if hx>hy: vf[d]=1
   for j in range(1,k): 
      W=LT_Wi(k,2*j,"Eq")
      U=random.sample(W,len(W)//2)
      for v in U: vf[LT_bin_to_int(v)]=1
   return BooleanFunction(list(vf))
   
   
   
   
   
   
