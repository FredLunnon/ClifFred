
################################################################################

#   "ClifFred" Python source for real degenerate Clifford algebras  Cl(p,q,r) . 
# Version 1.2; date: 13/02/18; author: Fred Lunnon <Fred.Lunnon@gmail.com> 
# In command window execute: 
#   cd /Users/fred/fred/python; python -i GA_multor.py 

from copy import *; 
import timeit; 

CFP_version = 1.2;  # update !! 

# Generate combination of  m  from  n  as array of  n  bits, with three sentinels; 
#   all subscripts from 0 ; terminate when false; universal ordering. 
#   Cases empty or  m = 0  fudged! 
def init_comb (n, m) :  # local j; 
  if n < 0 or m < 0 or m > n : comb = [];  # empty 
  else : comb = [ 0 for j in range(0, n-m) ] + [ 1 for j in range(0, m) ] + [1, 0, 0 ]; # end if 
  return comb; # end def 

def next_comb (comb) :  # local i,j,k,n3; 
  n3 = len(comb); 
  if n3 == 0 : return False;  # empty 
  i = 0; 
  while comb[i] == 0 : i = i+1; # end while 
  j = i; 
  while comb[j] == 1 : j = j+1; # end while 
  comb[j] = 1; 
  for k in range(0, j-i-1) : comb[k] = 1; # end for 
  for k in range(j-i-1, j) : comb[k] = 0; # end for 
  return j+1 <> n3 and j+3 <> n3; # end def 

# Binomial coefficient from formula for integer n,m : 
#   time = number of factors, or about O(2n); 
#   running answer is exact integer and not avoidably large; 
def binom (n, m) :  # local j, pro; 
  if m < 0 : return 0; 
  elif n >= 0 and m+m > n : return binom(n, n-m); 
  else : 
    pro = 1; 
    for j in range (0, m) : 
      pro = pro*(n-j)/(j+1); # end for 
    return pro; # end if end def 

class ClifFred : 
  
  def __init__(GA, siglis0 = [], siglis1 = []) :  
    # local j,k; # global eps,J; 
    GA.init(siglis0, siglis1);  # table set-up 
    GA.J = GA.bld([1]);  # quasi-pseudar 
    k = 0; 
    for j in range(0, GA.dim) : 
      if siglis0[j] <> 0 : 
        GA.J = GA.mul(GA.J, GA.gen(j+1)); 
        k = k+1; # end if end for 
    GA.J = GA.gra(GA.J, k); 
    GA.eps = 0.1**(16 - GA.dim); 
  # end def 
  
  # Sign of product monomial; flag false for Cl(n) case  
  def prosig (GA, moni, monj, flag = True) : 
    # local monh,moni,monj,sigh,sigi,t; 
    # global pow,dim,gensig,gramon; 
    sigh = +1;  # total sign of monh 
    sigi = +1 if GA.gramon[moni]%2 == 0 else -1;  #  t-th generator sign 
    for t in range(0, GA.dim) : 
      if moni & GA.pow[t] <> 0 : sigi = -sigi; # end if 
      if monj & GA.pow[t] <> 0 : 
        sigh = sigh*sigi; 
        if moni & GA.pow[t] <> 0 and flag : 
          sigh = sigh*GA.gensig[t]; # end if end if end for 
    return sigh; # end def 
  
  # Instantiate Clifford algebra, signature list comprising  +1,-1,0's 
  def init (GA, siglis0 = [], siglis1 = []) : 
    # local h,i,j,k,l,m,n,monh,moni,monj,t,sig,sig1; 
    # global dim,len,pow,bin,gensig,genlis,sigsig,mongo,gramon,offmon,siggo,doggo,leng3,ctrg3,ilisg3,jlisg3,hlisg3,slisg3; 
    GA.dim = len(siglis0);  # vector dimension 
    GA.bin = [ binom(GA.dim,k) for k in range(0, GA.dim+1) ];  # grade lengths 
    GA.pow = [ 2**k for k in range(0, GA.dim+1) ];  # binary powers 
    GA.len = GA.pow[GA.dim];  # multor dimension 
    
    # Generator signatures and their signs, extended optionally 
    GA.gensig = copy(siglis0); 
    GA.sigsig = siglis1 + [ siglis0[i] if siglis0[i] <> 0 else +1 if i >= len(siglis1) else siglis1[i] for i in range(len(siglis1), len(siglis0)) ]; 
    # Unit vectors  e1, .., en ?? 
    GA.genlis = [ [ [ ] ] + [ [ 1 if j == k else 0 for j in range(0, GA.dim) ] ] + [ [ ] for l in range(2, GA.dim+1) ] for k in range(0, GA.dim) ]; 

    # Monomial array and inverse, indexed by grade/offset 
    GA.mongo = [ [ None for i in range(0, GA.bin[k]) ] for k in range(0, GA.dim+1) ]; 
    GA.gramon = [ None for j in range(0, GA.len) ]; 
    GA.offmon = [ None for j in range(0, GA.len) ]; 
    for k in range(0, GA.dim+1) :  # grade  k  binary encoded monomials 
      i = 0;  # offset in grator component list 
      lis = init_comb(GA.dim, k);  # monomial qua binary list 
      while next_comb(lis) : 
        mon = 0;  # monomial qua binary integer 
        for t in range(0, GA.dim) : 
          mon = mon + lis[t]*GA.pow[t]; # end for 
        GA.mongo[k][i] = mon; GA.offmon[mon] = i; GA.gramon[mon] = k;  
        i = i+1; # end while end for 
    
    # Monomial signature and dual sign arrays, indexed by grade/offset 
    GA.siggo = [ [ None for i in range(0, GA.bin[k]) ] for k in range(0, GA.dim+1) ]; 
    GA.doggo = [ [ None for i in range(0, GA.bin[k]) ] for k in range(0, GA.dim+1) ]; 
    monh = GA.mongo[GA.dim][0];  # pseudar (spherical) 
    for k in range(0, GA.dim+1) : 
      for i in range(0, GA.bin[k]) : 
        moni = GA.mongo[k][i];  sig = 1; sig1 = 1; 
        for t in range(0, GA.dim) : 
          if moni & GA.pow[t] <> 0 : 
            sig = sig*GA.gensig[t]; sig1 = sig1*GA.sigsig[t]; # end if 
        GA.siggo[k][i] = sig; 
        GA.doggo[k][i] = sig1*GA.prosig(moni, monh, False); # end for end for 
    
    # Dry run: length and pter arrays for lists, indexed by grade triplet  k,l,m  
    GA.leng3 = [ [ [ 0 for k in range(0, GA.dim+1) ] for l in range(0, GA.dim+1) ] for m in range(0, GA.dim+1) ]; 
    GA.ctrg3 = [ [ [ 0 for k in range(0, GA.dim+1) ] for l in range(0, GA.dim+1) ] for m in range(0, GA.dim+1) ]; 
    for l in range(0, GA.dim+1) : 
      for k in range(0, GA.dim+1) :  
        for i in range(0, GA.bin[k]) : 
          moni = GA.mongo[k][i]; 
          for j in range(0, GA.bin[l]) : 
            monj = GA.mongo[l][j]; 
            # product monomial = +/- XOR(factor monomials) 
            m = GA.gramon[moni ^ monj]; 
            GA.leng3[m][l][k] = GA.leng3[m][l][k] + 1; 
            # end for end for end for end for end for 
    
    # Wet run: list of offset triplet and sign  i,j,h,s  arrays, indexed by  k,l,m  
    GA.ilisg3 = [ [ [ [ None for t in range(0, GA.leng3[m][l][k]) ] for k in range(0, GA.dim+1) ] for l in range(0, GA.dim+1) ] for m in range(0, GA.dim+1) ]; 
    GA.jlisg3 = [ [ [ [ None for t in range(0, GA.leng3[m][l][k]) ] for k in range(0, GA.dim+1) ] for l in range(0, GA.dim+1) ] for m in range(0, GA.dim+1) ]; 
    GA.hlisg3 = [ [ [ [ None for t in range(0, GA.leng3[m][l][k]) ] for k in range(0, GA.dim+1) ] for l in range(0, GA.dim+1) ] for m in range(0, GA.dim+1) ]; 
    GA.slisg3 = [ [ [ [ None for t in range(0, GA.leng3[m][l][k]) ] for k in range(0, GA.dim+1) ] for l in range(0, GA.dim+1) ] for m in range(0, GA.dim+1) ]; 
    for k in range(0, GA.dim+1) : 
      for l in range(0, GA.dim+1) :  
        for i in range(0, GA.bin[k]) : 
          moni = GA.mongo[k][i]; 
          for j in range(0, GA.bin[l]) : 
            monj = GA.mongo[l][j]; 
            # product monomial = +/- XOR(factor monomials) 
            monh = moni ^ monj; m = GA.gramon[monh]; h = GA.offmon[monh]; 
            t = GA.ctrg3[m][l][k]; 
            GA.ilisg3[m][l][k][t] = i; 
            GA.jlisg3[m][l][k][t] = j; 
            GA.hlisg3[m][l][k][t] = h; 
            GA.slisg3[m][l][k][t] = GA.prosig(moni, monj); 
            GA.ctrg3[m][l][k] = t + 1; 
            # end for end for end for end for end for 
    
    return None; # end def 
  
  
  # Extract  k-grator  <X>_k  of multor  X  qua list, expanding empty 
  def lis (GA, X, k = 0) :  # local h; # global dim,bin; 
    return [] if 0 > k or k > GA.dim else copy(X[k]) if X[k] <> 0 else [ 0 for h in range(0, GA.bin[m]) ]; # end def 
  
  # Grade  <X>_k  of multor  X  
  def gra (GA, X, k = 0) : 
    # local i,k,Z; global dim,bin; 
    Z = [ [ ] for i in range(0, GA.dim+1) ]; 
    if 0 <= k and k <= GA.dim : 
      Z[k] = copy(X[k]); # end if 
    return Z; # end def 
  
  # Grades of multor  X  from  lo  to  hi 
  def grades (GA, X, lo = None, hi = None) : 
    # local i,k,Z; global dim,bin; 
    if lo is None : lo = 0; hi = GA.dim;  # grade options 
    elif hi is None : hi = lo; # end if 
    Z = [ [ ] for i in range(0, GA.dim+1) ]; 
    for k in range(max(lo,0), min(hi,GA.dim)+1) : 
      Z[k] = copy(X[k]); # end if 
    return Z; # end def 
  
  # Even grades of multor  X  from  lo  to  hi 
  def even (GA, X, lo = None, hi = None) : 
    # local i,k,Z; global dim,bin; 
    if lo is None : lo = 0; hi = GA.dim;  # grade options 
    elif hi is None : hi = lo; # end if 
    Z = [ [ ] for i in range(0, GA.dim+1) ]; 
    for k in range(max(lo,0), min(hi,GA.dim)+1) : 
      if k % 2 == 0 : 
        Z[k] = copy(X[k]); # end if 
    return Z; # end def 
  
  # Odd grades of multor  X  from  lo  to  hi 
  def odd (GA, X, lo = None, hi = None) : 
    # local i,k,Z; global dim,bin; 
    if lo is None : lo = 0; hi = GA.dim;  # grade options 
    elif hi is None : hi = lo; # end if 
    Z = [ [ ] for i in range(0, GA.dim+1) ]; 
    for k in range(max(lo,0), min(hi,GA.dim)+1) : 
      if k % 2 <> 0 : 
        Z[k] = copy(X[k]); # end if 
    return Z; # end def 
  
  # Unit vector multor: return  Z = e_k , 1 <= k <= dim ; 
  def gen (GA, k) : 
    # local l,j; global dim,bin,genlis; 
    return copy(GA.genlis[k-1]); # end def 
  
  # Create multor from k-grator list  <X>_k ; list length check 
  def bld (GA, Xk_lis = [], k = 0) : 
    # local i,Z; global dim,bin; 
    assert 0 <= k and k <= GA.dim and (len(Xk_lis) == GA.bin[k] or len(Xk_lis) == 0), " grade or length wrong ";  # , k, GA.bin[k], len(Xk_lis)  o/p ?? 
    Z = [ [ ] for l in range(0, GA.dim+1) ]; 
    Z[k] = copy(Xk_lis); 
    return Z;  # end def 
  
  # Magnitude  ||X||  of multor  X 
  def mag2 (GA, X) : 
    # local i,k,s; global dim,bin,siggo; 
    s = 0; 
    for k in range(0, GA.dim+1) : 
      if X[k] <> [ ] : 
        for i in range(0, GA.bin[k]) :  # move up ?? 
          s = s + GA.siggo[k][i]*X[k][i]**2 
      # end if end for 
    return s; # end def 
  
  # Normalise multor  X , fails if isotropic : 
  #   requires square root, should test for small magnitude? 
  # def sta (GA, X) : 
  #   return GA.mul(X, GA.bld([1.0/sqrt(abs(GA.mag2(X)))])) # end def 
  
  # Reversion (main anti-automorphism)  X+  of multor  X 
  def rev (GA, X, lo = None, hi = None) : 
    # local i,k,sig,Z; global dim,bin; 
    if lo is None : lo = 0; hi = GA.dim;  # grade options 
    elif hi is None : hi = lo; # end if 
    Z = [ [ ] for k in range(0, GA.dim+1) ]; 
    for k in range(max(lo,0), min(hi,GA.dim)+1) : 
      sig = +1 if k%4 < 2 else -1; 
      if X[k] <> [ ] : 
        Z[k] = [ sig*X[k][i] for i in range(0, GA.bin[k]) ]; 
      # end if end for 
    return Z; # end def 
  
  # Parity transform (main automorphism)  X*  of multor  X 
  def par (GA, X, lo = None, hi = None) : 
    # local i,k,sig,Z; global dim,bin; 
    if lo is None : lo = 0; hi = GA.dim;  # grade options 
    elif hi is None : hi = lo; # end if 
    Z = [ [ ] for k in range(0, GA.dim+1) ]; 
    for k in range(max(lo,0), min(hi,GA.dim)+1) : 
      sig = +1 if k%2 < 1 else -1; 
      if X[k] <> [ ] : 
        if sig > 0 : Z[k] = copy(X[k]); 
        else : Z[k] = [ -X[k][i] for i in range(0, GA.bin[k]) ]; 
      # end if end for 
    return Z; # end def 
  
  # Dual  Z  of multor  X  (combinatorial: null signatures signed!) 
  def dual (GA, X, lo = None, hi = None) : 
    # local i,j,k,Z; global dim,bin,siggo,doggo; 
    if lo is None : lo = 0; hi = GA.dim;  # grade options 
    elif hi is None : hi = lo; # end if 
    Z = [ [ ] for l in range(max(GA.dim-hi,0), min(GA.dim-lo,GA.dim)+1) ]; 
    for k in range(0, GA.dim+1) : 
      if X[k] <> [ ] : 
        j = GA.bin[k]-1;  #  j = GA.bin[m] 
        Z[GA.dim-k] = [ GA.doggo[k][j-i]*X[k][j-i] for i in range(0, j+1) ]; 
    # end if end for 
    return Z;  # end def 
  
  # Add multors: return  Z = X + Y 
  def add (GA, X, Y, lo = None, hi = None) : 
    # local i,k,Z; global dim,bin; 
    if lo is None : lo = 0; hi = GA.dim;  # grade options 
    elif hi is None : hi = lo; # end if 
    Z = [ [ ] for k in range(0, GA.dim+1) ]; 
    for k in range(max(lo,0), min(hi,GA.dim)+1) : 
      if X[k] == [ ] : Z[k] = copy(Y[k]); 
      elif Y[k] == [ ] : Z[k] = copy(X[k]); 
      else : Z[k] = [ X[k][i] + Y[k][i] for i in range(0, GA.bin[k]) ]; 
      # end if end for 
    return Z; # end def 
  
  # Subtract multors: return  Z = X - Y 
  def sub (GA, X, Y, lo = None, hi = None) : 
    # local i,k,Z; global dim,bin; 
    if lo is None : lo = 0; hi = GA.dim;  # grade options 
    elif hi is None : hi = lo; # end if 
    Z = [ [ ] for k in range(0, GA.dim+1) ]; 
    for k in range(max(lo,0), min(hi,GA.dim)+1) : 
      if X[k] == [ ] : 
        Z[k] = [ ] if (Y[k] == [ ]) else [ -Y[k][i] for i in range(0, GA.bin[k]) ]; 
      elif Y[k] == [ ] : Z[k] = copy(X[k]); 
      else : Z[k] = [ X[k][i] - Y[k][i] for i in range(0, GA.bin[k]) ]; 
      # end if end for 
    return Z; # end def 
  
  # Multiply multors: return  Z = X Y , grades  lo  to  hi 
  def mul (GA, X, Y, lo = None, hi = None) : 
    # global dim,bin,leng3,ilisg3,jlisg3,hlisg3,slisg3; 
    # local i,j,k,l,m,h,s,t,Z,Zmlk,hlis,ilis,jlis,slis; 
    if lo is None : lo = 0; hi = GA.dim;  # grade options 
    elif hi is None : hi = lo; # end if 
    Z = [ [ ] for m in range(0, GA.dim+1) ]; 
    for m in range(max(lo,0), min(hi,GA.dim)+1) : 
      for l in range(0, GA.dim+1) :  # favour  Y  shorter than  X  
        if Y[l] <> [ ] : 
          for k in range(abs(m-l), min(m+l, GA.dim+GA.dim-m-l)+1, 2) : 
            if X[k] <> [ ] and GA.leng3[m][l][k] > 0 : 
              ilis = GA.ilisg3[m][l][k]; 
              jlis = GA.jlisg3[m][l][k]; 
              hlis = GA.hlisg3[m][l][k]; 
              slis = GA.slisg3[m][l][k]; 
              if Z[m] == [ ] : Z[m] = [ 0 for h in range(0, GA.bin[m]) ]; 
              for t in range(0, GA.leng3[m][l][k]) : 
                Z[m][hlis[t]] = Z[m][hlis[t]] + slis[t]*X[k][ilis[t]]*Y[l][jlis[t]]; 
              # end for end if end for end if end for end for 
    return Z; # end def 
  
  # Dual multiply: return  Z = (X~ Y~)~ , grades  lo  to  hi 
  def lum (GA, X, Y, lo = None, hi = None) : 
    if lo is None : lo = 0; hi = GA.dim;  # grade options 
    elif hi is None : hi = lo; # end if 
    return GA.dual(GA.mul(GA.dual(X), GA.dual(Y), GA.dim-hi, GA.dim-lo));  # end def 
  
  # Naive left contract: return  Z = X _| Y = \sum_k,l < <X>_k <Y>_l >_l-k  
  def lecon0 (GA, X, Y) : 
    # local k,l,Z; 
    Z = GA.bld([0]); 
    for l in range(0, GA.dim+1) : 
      for k in range(0, l+1) : 
        Z = GA.add(Z, GA.mul(GA.gra(X, k), GA.gra(Y, l), l-k)); 
    # end for end for 
    return Z; # end def 
  
  # Naive right contract: return  Z = X |_ Y = \sum_k,l < <X>_k <Y>_l >_k-l  
  def ricon0 (GA, X, Y) : 
    # local k,l,Z; 
    Z = GA.bld([0]); 
    for l in range(0, GA.dim+1) : 
      for k in range(l, GA.dim+1) : 
        Z = GA.add(Z, GA.mul(GA.gra(X, k), GA.gra(Y, l), k-l)); 
    # end for end for 
    return Z; # end def 
  
  # Naive outer product: return  Z = X ^ Y = \sum_k,l < <X>_k <Y>_l >_k+l  
  def wedge0 (GA, X, Y) : 
    # local k,l,Z; 
    Z = GA.bld([0]); 
    for l in range(0, GA.dim+1) : 
      for k in range(0, GA.dim-l+1) : 
        Z = GA.add(Z, GA.mul(GA.gra(X, k), GA.gra(Y, l), k+l)); 
    # end for end for 
    return Z; # end def 
  
  # Naive Dorst inner product: return  Z = X o Y = \sum_k,l < <X>_k <Y>_l >_|k-l|  
  def fat0 (GA, X, Y) : 
    # local k,l,Z; 
    Z = GA.bld(); 
    for l in range(0, GA.dim+1) : 
      for k in range(0, GA.dim+1) : 
        Z = GA.add(Z, GA.mul(GA.gra(X, k), GA.gra(Y, l), abs(k-l))); 
    # end for end for 
    return Z; # end def 
  
  # Naive Hestenes inner product: return  Z = X o Y = \sum_k,l < <X>_k <Y>_l >_|k-l|  
  def dot0 (GA, X, Y) : 
    # local k,l,Z; 
    Z = GA.bld(); 
    for l in range(1, GA.dim+1) : 
      for k in range(1, GA.dim+1) : 
        Z = GA.add(Z, GA.mul(GA.gra(X, k), GA.gra(Y, l), abs(k-l))); 
    # end for end for 
    return Z; # end def 
  
  # Naive scalar product: return  Z = X * Y = \sum_k,l < <X>_k <Y>_l >_0  
  def scalar0 (GA, X, Y) : 
    # local k,l,Z; 
    Z = GA.bld(); 
    for l in range(0, GA.dim+1) : 
      for k in range(0, GA.dim+1) : 
        Z = GA.add(Z, GA.mul(GA.gra(X, k), GA.gra(Y, l), 0)); 
    # end for end for 
    return Z; # end def 
  
  # Left contract multors: return  Z = X _| Y , grades  lo  to  hi 
  #   \sum_k,l < <X>_k <Y>_l >_l-k , so  l >= m  &  k = l-m 
  def lecon (GA, X, Y, lo = None, hi = None) : 
    # global dim,bin,leng3,ilisg3,jlisg3,hlisg3,slisg3; 
    # local i,j,k,l,m,h,s,t,Z,Zmlk,hlis,ilis,jlis,slis; 
    if lo is None : lo = 0; hi = GA.dim;  # grade options 
    elif hi is None : hi = lo; # end if 
    Z = [ [ ] for m in range(0, GA.dim+1) ]; 
    for m in range(max(lo,0), min(hi,GA.dim)+1) : 
      for l in range(m, GA.dim+1) :  # favour  Y  shorter than  X  
        if Y[l] <> [ ] : 
          for k in [l-m] :  #  l >= m  &  k = l-m 
            if X[k] <> [ ] and GA.leng3[m][l][k] > 0 : 
              ilis = GA.ilisg3[m][l][k]; 
              jlis = GA.jlisg3[m][l][k]; 
              hlis = GA.hlisg3[m][l][k]; 
              slis = GA.slisg3[m][l][k]; 
              if Z[m] == [ ] : Z[m] = [ 0 for h in range(0, GA.bin[m]) ]; 
              for t in range(0, GA.leng3[m][l][k]) : 
                Z[m][hlis[t]] = Z[m][hlis[t]] + slis[t]*X[k][ilis[t]]*Y[l][jlis[t]]; 
              # end for end if end for end if end for end for 
    return Z; # end def 
  
  # Right contract multors: return  Z = X _| Y , grades  lo  to  hi 
  #   \sum_k,l < <X>_k <Y>_l >_k-l  so  l <= dim-m  &  k = l+m 
  def ricon (GA, X, Y, lo = None, hi = None) : 
    # global dim,bin,leng3,ilisg3,jlisg3,hlisg3,slisg3; 
    # local i,j,k,l,m,h,s,t,Z,Zmlk,hlis,ilis,jlis,slis; 
    if lo is None : lo = 0; hi = GA.dim;  # grade options 
    elif hi is None : hi = lo; # end if 
    Z = [ [ ] for m in range(0, GA.dim+1) ]; 
    for m in range(max(lo,0), min(hi,GA.dim)+1) : 
      for l in range(0, GA.dim-m+1) :  # l <= dim-m  
        if Y[l] <> [ ] : 
          for k in [m+l] :  #  k = m+l 
            if X[k] <> [ ] and GA.leng3[m][l][k] > 0 : 
              ilis = GA.ilisg3[m][l][k]; 
              jlis = GA.jlisg3[m][l][k]; 
              hlis = GA.hlisg3[m][l][k]; 
              slis = GA.slisg3[m][l][k]; 
              if Z[m] == [ ] : Z[m] = [ 0 for h in range(0, GA.bin[m]) ]; 
              for t in range(0, GA.leng3[m][l][k]) : 
                Z[m][hlis[t]] = Z[m][hlis[t]] + slis[t]*X[k][ilis[t]]*Y[l][jlis[t]]; 
              # end for end if end for end if end for end for 
    return Z; # end def 
  
  # Outer product of multors: return  Z = X ^ Y , grades  lo  to  hi 
  #   \sum_k,l < <X>_k <Y>_l >_k+l  so  l <= dim-k  &  k = m-l 
  def wedge (GA, X, Y, lo = None, hi = None) : 
    # global dim,bin,leng3,ilisg3,jlisg3,hlisg3,slisg3; 
    # local i,j,k,l,m,h,s,t,Z,Zmlk,hlis,ilis,jlis,slis; 
    if lo is None : lo = 0; hi = GA.dim;  # grade options 
    elif hi is None : hi = lo; # end if 
    Z = [ [ ] for m in range(0, GA.dim+1) ]; 
    for m in range(max(lo,0), min(hi,GA.dim)+1) : 
      for l in range(0, m+1) :  # l <= m 
        if Y[l] <> [ ] : 
          for k in [m-l] :  #  k = m-l 
            if X[k] <> [ ] and GA.leng3[m][l][k] > 0 : 
              ilis = GA.ilisg3[m][l][k]; 
              jlis = GA.jlisg3[m][l][k]; 
              hlis = GA.hlisg3[m][l][k]; 
              slis = GA.slisg3[m][l][k]; 
              if Z[m] == [ ] : Z[m] = [ 0 for h in range(0, GA.bin[m]) ]; 
              for t in range(0, GA.leng3[m][l][k]) : 
                Z[m][hlis[t]] = Z[m][hlis[t]] + slis[t]*X[k][ilis[t]]*Y[l][jlis[t]]; 
              # end for end if end for end if end for end for 
    return Z; # end def 
  
  # Dual outer product: return  Z = X v Y = (X~ ^ Y~)~ 
  def vee (GA, X, Y, lo = None, hi = None) : 
    if lo is None : lo = 0; hi = GA.dim;  # grade options 
    elif hi is None : hi = lo; # end if 
    return GA.dual(GA.wedge(GA.dual(X), GA.dual(Y), GA.dim-hi, GA.dim-lo));  # end def 
  
  # Dorst inner product multors: return  Z = X o Y , grades  lo  to  hi 
  #   \sum_k,l < <X>_k <Y>_l >_|k-l|  so  l <= dim-m  &  k = l+m  or  l >= m  &  k = l-m ;  once only when  m = 0 ! 
  def fat (GA, X, Y, lo = None, hi = None) : 
    # global dim,bin,leng3,ilisg3,jlisg3,hlisg3,slisg3; 
    # local i,j,k,l,m,h,s,t,Z,Zmlk,hlis,ilis,jlis,slis; 
    if lo is None : lo = 0; hi = GA.dim;  # grade options 
    elif hi is None : hi = lo; # end if 
    Z = [ [ ] for m in range(0, GA.dim+1) ]; 
    for m in range(max(lo,0), min(hi,GA.dim)+1) : 
      for l in range(0, GA.dim+1) :  #  l <= dim-m ; l >= m  
        if Y[l] <> [ ] : 
          for k in ([l+m] if l <= GA.dim-m else []) + ([l-m] if l >= m and m > 0 else []) :  #  k = {l+m} U {l-m} 
            if X[k] <> [ ] and GA.leng3[m][l][k] > 0 : 
              ilis = GA.ilisg3[m][l][k]; 
              jlis = GA.jlisg3[m][l][k]; 
              hlis = GA.hlisg3[m][l][k]; 
              slis = GA.slisg3[m][l][k]; 
              if Z[m] == [ ] : Z[m] = [ 0 for h in range(0, GA.bin[m]) ]; 
              for t in range(0, GA.leng3[m][l][k]) : 
                Z[m][hlis[t]] = Z[m][hlis[t]] + slis[t]*X[k][ilis[t]]*Y[l][jlis[t]]; 
              # end for end if end for end if end for end for 
    return Z; # end def 
  
  # Hestenes inner product multors: return  Z = X . Y , grades  lo  to  hi 
  #   \sum_k>0,l>0 < <X>_k <Y>_l >_|k-l|  so  l <= dim-m  &  k = l+m  or  l >= m  &  k = l-m ;  m > 0 ! 
  def dot (GA, X, Y, lo = None, hi = None) : 
    # global dim,bin,leng3,ilisg3,jlisg3,hlisg3,slisg3; 
    # local i,j,k,l,m,h,s,t,Z,Zmlk,hlis,ilis,jlis,slis; 
    if lo is None : lo = 0; hi = GA.dim;  # grade options 
    elif hi is None : hi = lo; # end if 
    Z = [ [ ] for m in range(0, GA.dim+1) ]; 
    for m in range(max(lo,0), min(hi,GA.dim)+1) : 
      for l in range(1, GA.dim+1) :  #  l <= dim-m ; l >= m  
        if Y[l] <> [ ] : 
          for k in ([l+m] if l <= GA.dim-m else []) + ([l-m] if l > m and m > 0 else []) :  #  k = {l+m} U {l-m} 
            if X[k] <> [ ] and GA.leng3[m][l][k] > 0 : 
              ilis = GA.ilisg3[m][l][k]; 
              jlis = GA.jlisg3[m][l][k]; 
              hlis = GA.hlisg3[m][l][k]; 
              slis = GA.slisg3[m][l][k]; 
              if Z[m] == [ ] : Z[m] = [ 0 for h in range(0, GA.bin[m]) ]; 
              for t in range(0, GA.leng3[m][l][k]) : 
                Z[m][hlis[t]] = Z[m][hlis[t]] + slis[t]*X[k][ilis[t]]*Y[l][jlis[t]]; 
              # end for end if end for end if end for end for 
    return Z; # end def 
  
  # Scalar product multors: return  Z = X * Y 
  #   \sum_k,l < <X>_k <Y>_l >_0 , so  m = 0  &  k = l 
  def scalar (GA, X, Y) : 
    # global dim,bin,leng3,ilisg3,jlisg3,hlisg3,slisg3; 
    # local i,j,k,l,m,h,s,t,Z,Zmlk,hlis,ilis,jlis,slis; 
    Z = [ [ ] for m in range(0, GA.dim+1) ]; 
    for m in [0] : 
      for l in range(0, GA.dim+1) :  # favour  Y  shorter than  X  
        if Y[l] <> [ ] : 
          for k in [l] :  #  l >= m  &  k = l-m 
            if X[k] <> [ ] and GA.leng3[m][l][k] > 0 : 
              ilis = GA.ilisg3[m][l][k]; 
              jlis = GA.jlisg3[m][l][k]; 
              hlis = GA.hlisg3[m][l][k]; 
              slis = GA.slisg3[m][l][k]; 
              if Z[m] == [ ] : Z[m] = [ 0 for h in range(0, GA.bin[m]) ]; 
              for t in range(0, GA.leng3[m][l][k]) : 
                Z[m][hlis[t]] = Z[m][hlis[t]] + slis[t]*X[k][ilis[t]]*Y[l][jlis[t]]; 
              # end for end if end for end if end for end for 
    return Z; # end def 
  
  # Conjugation transform of multor: return  Z = (Y+) X Y , grades  lo  to  hi  
  def form (GA, X, Y, lo = None, hi = None) : 
    if lo is None : lo = 0; hi = GA.dim;  # grade options 
    elif hi is None : hi = lo; # end if 
    return GA.mul(GA.mul(GA.rev(Y), X), Y, lo, hi);  # end def 
  
  # Sum list of multors 
  def addlis (GA, Ylis) :  # local Z,i; 
    Z = GA.bld([0]); 
    for i in range(0, len(Ylis)) : Z = GA.add(Z, Ylis[i]); # end for 
    return Z; # end def 
  
  # Product list of multors 
  def mullis (GA, Ylis) :  # local Z,i; 
    Z = GA.bld([1]); 
    for i in range(0, len(Ylis)) : Z = GA.mul(Z, Ylis[i]); # end for 
    return Z; # end def 
  
  # Retrieve component of monomial in  X 
  def getcom (GA, X, mon = 0) :  # local k, Xk; 
    assert 0 < mon or mon > GA.len, " illegal monomial ";  
    k = GA.gramon[mon]; 
    return X[k][GA.offmon[mon]] if k >= 0 and X[k] <> [] else 0; # end def 
  
  # Set component of monomial in  X 
  def putcom (GA, X, mon = 0, val = 0) :  # local k; 
    assert 0 < mon or mon > GA.len, " illegal monomial ";  
    k = GA.gramon[mon];  # now build if empty 
    if X[k] == [] : X[k] = [ 0 for i in range(0, GA.bin[k]) ]; # end if 
    X[k][GA.offmon[mon]] = val; 
    return None; # end def 
  
  # Grade  k  and offset  t  of monomial 
  def graset (GA, mon) : 
    assert 0 < mon or mon > GA.len, " illegal monomial ";  
    return GA.gramon[mon], GA.offmon[mon]; # end def 
  
  # Monomial at grade  k  and offset  t 
  def monom (GA, k, t) : 
    return GA.mongo[k][t] if 0 <= k and k <= GA.dim and 0 <= t and t < GA.bin[k]else 0; # end def 
  
  
  # Distance from multor  X  to zero 
  def zero (GA, X) : 
    # local k, t, err; global GA.bin, GA.dim; 
    err = 0.0;  # error 
    for k in range(0, GA.dim+1) : 
      if X[k] <> [] : 
        for t in range(0, GA.bin[k]) : 
          err = err + abs(X[k][t]); # end for end if end for 
    err = err/GA.len; 
    return err; # end def 
  
  # Distance from multors  X, Y  to proportionality 
  def prop (GA, X, Y) : 
    # local k, t, err, Xcs, Ycs; global GA.bin, GA.dim; 
    Xcs = 0.0;  Ycs = 0.0;  # component sums 
    for k in range(0, GA.dim+1) : 
      for t in range(0, GA.bin[k]) : 
        if X[k] <> [] : Xcs = Xcs + abs(X[k][t]); # end if 
        if Y[k] <> [] : Ycs = Ycs + abs(Y[k][t]); # end if end for end for 
    err = 0.0; 
    for k in range(0, GA.dim+1) : 
      for t in range(0, GA.bin[k]) : 
        if X[k] <> [] and Y[k] <> [] : err = err + abs(Xcs*Y[k][t] - Ycs*X[k][t]); 
        elif X[k] <> [] : err = err + abs(Ycs*X[k][t]); 
        elif Y[k] <> [] : err = err + abs(Xcs*Y[k][t]); # end if end for end for 
    err = err/(2.0*Xcs*Ycs) if (Xcs+Ycs)/2.0 > GA.eps else 1.0; 
    return err; # end def 
  
  # Is multor  X  near zero? 
  def is_zero (GA, X) : 
    # global GA.eps; 
    return GA.zero(X) < GA.eps; # end def 
  
  # Are multors  X, Y  nearly proportional? 
  def is_prop (GA, X, Y) : 
    # global GA.eps; 
    return GA.prop(X, Y) < GA.eps; # end def 
  
  # end class 

multor = ClifFred;  # class variable for wrapper switching 

if False : # if True :  # visual check tables 
  GA = ClifFred([+1,+1]); 
  print GA.dim; print GA.len; print GA.pow; print GA.bin; 
  print GA.gensig; print GA.sigsig; print GA.mongo; print GA.gramon; print GA.offmon; 
  print GA.siggo; print GA.doggo; print GA.leng3; print GA.ctrg3; 
  print GA.ilisg3; print GA.jlisg3; print GA.hlisg3; print GA.slisg3; 
# end if 

# Demo integer test of algebraic identities:  Cl(6) , Cl(3,3), Cl(2,2,2) ; 
print "ClifFred testing version ", CFP_version; 
secs = timeit.default_timer(); 
for (sigs, mags) in \
[ ( [+1,+1,+1,+1,+1,+1], [30517578125, 15625, 12064000] ), \
( [+1,+1,+1,-1,-1,-1], [307546875, 3375, -4193280] ), \
( [+1,+1,-1,-1,+0,+0], [530841600, 3600, 38400] ) ] : 
  GA = multor(sigs);  # initialise algebra 
  X = GA.mullis([ GA.add(GA.bld([2]), GA.mul(GA.gen(j+1), GA.gen(i+1))) \
  for i in range(0, GA.dim) for j in range(0, i) ]); X;  
  # X = (2 + e1 e2)(2 + e1 e3)(2 + e2 e3) ... 
  Y = GA.mullis([ GA.add(GA.bld([2]), GA.gen(i+1)) \
  for i in range(0, GA.dim) ]); Y;  # Y = (2 + e1)(2 + e2)(2 + e3) ... 
  Z = GA.mullis([ GA.add(GA.bld([2]), GA.mul(GA.bld([i+1]), GA.gen(i+1))) \
  for i in range(0, GA.dim) ]); Z;  # Z = (2 + e1)(2 + 2 e2)(2 + 3 e3) ... 
  
  assert GA.add(GA.even(Z), GA.odd(Z)) == Z, \
  "ClifFred:  even(Z) + odd(Z) = Z "; 
  assert GA.add(GA.grades(Z, 0, 3), GA.grades(Z, 4, 6)) == Z, \
  "ClifFred:  <Z>_{0,3} + <Z>_{4,6} = Z "; 
  assert [GA.mag2(X), GA.mag2(Y), GA.mag2(Z)] == mags, \
  "ClifFred:  ||X||, ||Y||, ||Z|| "; 
  assert GA.mag2(X) == GA.lis(GA.mul(GA.rev(X), X))[0], \
  "ClifFred:  ||X|| = < X+ X >_0 "; 
  assert GA.mul(GA.mul(X, Y), Z) == GA.mul(X, GA.mul(Y, Z)), \
  "ClifFred:  X (Y Z) = (X Y) Z "; 
  assert GA.mul(X, GA.add(Y, Z)) == GA.add(GA.mul(X, Y), GA.mul(X, Z)), \
  "ClifFred:  X (Y + Z) = X Y + X Z "; 
  assert GA.lecon(GA.wedge(X, Y), Z) == GA.lecon(X, GA.lecon(Y, Z)), \
  "ClifFred:  (X ^ Y) _| Z =  X _| (Y _| Z) "; 
  assert GA.add(GA.fat(X, Z), GA.scalar(X, Z)) == GA.add(GA.lecon(X, Z), GA.ricon(X, Z)), \
  "ClifFred:  (X o Z) + (X * Z) = (X _| Z) + (X |_ Z) "; 
  assert GA.add(GA.fat(X, Z), GA.mul(GA.gra(X), GA.gra(Z))) \
  == GA.add(GA.dot(X, Z), GA.add(GA.mul(GA.gra(X), Z), GA.mul(X, GA.gra(Z)))), \
  "ClifFred:  (X o Z) + <X>_0 <Z>_0 = (X . Z) + <X>_0 Z  + X <Z>_0 "; 
  assert GA.rev(GA.wedge(Z, Y)) == GA.wedge(GA.rev(Y), GA.rev(Z)), \
  "ClifFred:  (Z ^ Y)+ = (Y+) ^ (Z+) "; 
  assert GA.wedge(X, GA.wedge(Y, Z)) == GA.wedge(GA.wedge(X, Y), Z), \
  "ClifFred:  X ^ (Y ^ Z) = (X ^ Y) ^ Z "; 
  assert GA.vee(X, GA.vee(Y, Z)) == GA.vee(GA.vee(X, Y), Z), \
  "ClifFred:  X v (Y v Z) = (X v Y) v Z "; 
  assert GA.rev(GA.mul(Z, Y)) == GA.mul(GA.rev(Y), GA.rev(Z)), \
  "ClifFred:  (Z Y)+ = (Y+) (Z+) "; 
  assert GA.par(GA.mul(Z, Y)) == GA.mul(GA.par(Z), GA.par(Y)), \
  "ClifFred:  (Z Y)* = (Z*) (Y*) "; 
  assert GA.rev(GA.rev(Z)) == Z, \
  "ClifFred:  (Z+)+ = Z "; 
  assert GA.par(GA.par(Z)) == Z, \
  "ClifFred:  (Z*)* = Z "; 
  s = 1 - ( GA.dim - sum(GA.sigsig)//2 )%2*2;  # sign change via squared dual 
  assert GA.dual(GA.dual(Z)) == GA.mul(GA.bld([s]), Z), \
  "ClifFred:  (Z~)~ = +/- Z "; 

# end for 
secs = timeit.default_timer() - secs; 
print " elapsed time in secs ", secs;  # 0.13 sec 


# TODO --- 
# (versor) divide, pre-divide; rescale by 1/integer; normalise() with sqrt; 
#   Attach  GA.s0, GA.s1, GA.e[k]  --- ?? 
# GA.neg() subtract from GA.s0 = [ [ ] for k in range(0, GA.dim+1) ]; 
# Maybe sum, product, wedge, ... accept arbitrarily many args?  
#   init() etc :  h -> k -> l -> m -> n  ... 
#   class ClifFred (object) : https://docs.python.org/2/tutorial/classes.html sect. 9.5 inheritance. 
# Commit --- 
#   Wrapper version and time printed by short test 
#   is_zero(), is_prop() : Boolean versions 
#   GA.dual() now graded, as documentation 
#   lis() now expands empty to zeros 
#   genlis[] now restored 
#   bld() : now checks list length 
#   GA.vee() fixed & tested (Kurt Nalty) 
#   Renamed GA.dot() as GA.fat() ; Hestenes GA.dot() now 
#   Renamed ClifFred() as multor() ; switchable between wrappers 

# Polynomial I/O; default identifiers (SymPy) ? 
#   http://mattpap.github.io/scipy-2011-tutorial/html/basics.html 
# Overload arithmetic +,-,*  etc. ? 
#   http://www.programiz.com/python-programming/operator-overloading 

################################################################################

