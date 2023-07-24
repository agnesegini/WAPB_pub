def h16():
  """Return the hybrid Boolean function in 16 variables h_16 of Section 4.3 [GM22b]"""
  f=CMR(16)
  F=LM8()
  G=construction_2(8,F)
  n=16
  nG= rev_fun(G)
  nf= rev_fun(f)
  H=Ikn(0,n)*f+Ikn(1,n)*nf+Ikn(2,n)*nf
  for k in range(3,8):
    H+=Ikn(k,n)*nG
  for k in range(8,14):
    H+=Ikn(k,n)*G
  for k in range(14,17):
    H+=Ikn(k,n)*f
  return H
