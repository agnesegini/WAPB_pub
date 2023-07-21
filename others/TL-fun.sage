"""
Building WAPB functions from https://dl.acm.org/doi/abs/10.1007/s12095-019-00374-6 

To obtain a WPB function the construction as to be initialised with a sequence c=[1,0,0,...,0,-1] (e.g, c=[1,0,0,0,-1]).
See Remark Remark 1 of https://ia.cr/2023/110. 
"""


def LT_Wi(k,i,rel):
   """Returns  W^{rel}_(2k,i) for LT constructions, i.e. subset of the set
   of all vectors (x,y) in F2^{2k} with Hamming weight i and such that Hw(x)~rel~Hw(y).
   

   Args:
       k (int): 2k is the size of the vectors
       i (int): hamming weight
       rel (string): "Eq" for = and "Gr" for >.

   Returns:
       list:  the list of all elements of W^{rel}_(2k,i)
   """   
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
   """
   Returns the set W^> the union of W^{>}_(2k,i) for i=0,...,2k.

   Args:
       k (int): 2k is the size of the vectors

   Returns:
       set: W^>
   """   
   n=2*k 
   W=set(LT_Wi(k,0,"Gr"))
   for i in range(0,n+1):
      Wi=LT_Wi(k,i,"Gr")
      for w in Wi: W.add(w)
   return W
       
def LT_Ui(k,i,c):
   """Returns an arbitrary subset  Ui of W^{=}_(2k,i) of cardinality (|W^{=}_(2k,i)|+c)/2

   Args:
       k (int): 2k is the size of the vectors
       i (int): hamming weight
       c (int): c in {0,1,-1} is the cardinality corrections 

   Returns:
       set: Ui
   """   
   if mod(i,2)==1: return [] 
   else: 
      W=LT_W(k,i,"Eq")
      w=len(W)
      ui=(w+c)//2
      U=random.sample(W, ui)
      assert len(U)==ui
      return U

def LT_U(k,c):
   """Returns the union U of subsets  U_{2i} of W^{=}_(2k,2i) of cardinality (|W^{=}_(2k,2i)|+c[i])/2
      for k=0..k+1

   Args:
       k (int): 2k is the size of the vectors
       c (list): sequence of lenght k+1 over {−1, 0, 1} suck that c[i] = 0 if {k choose i} is event, and c[i] = +/-1 if odd. 

   Returns:
       set: U
   """   
   n=2*k 
   U=set(LT_Ui(k,0,c[0]))
   for i in range(1,k+1):
      U2i=LT_Ui(k,2*i,c[i])
      for w in U2i: U.add(w)
   return U

def LT_is_valid_seq(c,k):
   """Checks if c is a valid input sequence for LT_U(k,c).

   Args:
       c (list): sequence over {−1, 0, 1} of lenght k+1
       k (int): 2k is the size of the vectors

   Returns:
       bool: True iff valid 
   """   
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
   """Returns a boolean function function in n=2k variables from TL construction with a given sequence as input

   Args:
       c (list): sequence of lenght k+1 over {−1, 0, 1} suck that c[i] = 0 if {k choose i} is event, and c[i] = +/-1 if odd. 
       k (int): 2k is the size of the vectors

   Returns:
       BooleanFunction: function from TL family
   """   
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
   """
   Returns a WPB function in n variables sampled uniformly at random from TL family

   Args:
       n (int): number of variables

   Returns:
       BooleanFunction: TL-Boolean function in n variables that is WPB
   """   
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
   
   
   
   
   
   
