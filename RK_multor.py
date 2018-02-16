
################################################################################

#   <clifford> compatible Python wrapper for Robert Kern <clifford> package. 
# Version 0; date: 15/02/18; author: Fred Lunnon <Fred.Lunnon@gmail.com> 
# In command window execute: 				 --- PARTIAL ?? 
#   cd /Users/fred/fred/python; python -i RK_multor.py 
# Links and file paths both fail; scratch clifford-master/clifford copy! 

from copy import *; 
import timeit; 
from clifford import *; 
from clifford import _eps;  #  --- FAILS ?? 

class multor : 
  
  def __init__(GA, siglis0 = [], siglis1 = []) : 
    # local p, q, k, GAlayout, GAblades; 
    # global dim,len,pow,bin,binsum,gensig,sigsig,genlis,_eps,layout,bladesJ; 
    global _eps;  #  --- FAILS ?? 
    GA.dim = len(siglis0);  # vector dimension 
    GA.bin = [ comb(GA.dim,k) for k in range(0, GA.dim+1) ];  # grade lengths 
    GA.binsum = [ sum(GA.bin[0:k]) for k in range(0, GA.dim+2) ];  # grade offsets 
    GA.pow = [ 2**k for k in range(0, GA.dim+1) ];  # binary powers 
    GA.len = GA.pow[GA.dim];  # multor dimension 
    GA.eps = 0.1**(16 - GA.dim); _eps = GA.eps;  #  --- FAILS ?? 
    
    GA.gensig = copy(siglis0); GA.sigsig = copy(siglis0); 
    p = 0; q = 0;  # enforce consecutive signs in {+1,-1} 
    for k in range(0, GA.dim) : 
      if GA.gensig[k] == +1 and q == 0 : p = p+1; 
      elif GA.gensig[k] == -1 : q = q+1; 
      else : assert False, " unsupported signature "; # end if end for 
    GA.layout, GA.blades = Cl(p, q);  # RK initialiser 
    GA.genlis = [ GA.blades[GA.layout.names[k+1]] for k in range(0, GA.dim)]; 
    # disinter unit vectors  e1,..,en : is there another way? 
    # polynomial print names are  'e0', ...,  not  'e1', ... ?! 
    
    GA.J = GA.bld([1]);  # quasi-pseudar 
    k = 0; 
    for j in range(0, GA.dim) : 
      if siglis0[j] <> 0 : 
        GA.J = GA.mul(GA.J, GA.gen(j+1)); 
        k = k+1; # end if end for 
    GA.J = GA.gra(GA.J, k); 
  # end def 
  
  # Extract  k-grator  <X>_k  of multor  X  qua list 
  def lis (GA, X, k = 0) :  # global dim,bin; 
    return X.value[GA.binsum[k]:GA.binsum[k+1]] if 0 <= k and k <= GA.dim else [];  # end def 
  
  # Grade  <X>_k  of multor  X 
  def gra (GA, X, k = 0) : 
    # global dim,bin; 
    if 0 <= k and k <= GA.dim : Z = X.__call__(k); 
    else : Z = MultiVector(GA.layout); # end if 
    return Z; # end def 
  
  # Grades  <X>_k  of multor  X  for  lo <= k <= hi 
  def grades (GA, X, lo, hi) : 
    # local h,hk,hl,Z; global dim,binsum,layout; 
    hk = GA.binsum[max(lo, 0)]; hl = GA.binsum[min(hi, GA.dim)+1]; 
    Z = MultiVector(GA.layout); 
    Z.value[hk:hl] = X.value[hk:hl]; 
    return Z; # end def 
  
  # Even grades of multor  X  from  lo  to  hi 
  def even (GA, X, lo = None, hi = None) : 
    # local i,k,Z; global dim,bin; 
    assert False, " UNIMPLEMENTED ";  #  --- TEMP ?? 
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
    assert False, " UNIMPLEMENTED ";  #  --- TEMP ?? 
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
  
  # Create multor from k-grator list  <X>_k ; list length check; RK omitted! 
  def bld (GA, Xk_lis = [], k = 0) : 
    # local h,h0,Z; global dim,bin,binsum,layout; 
    assert 0 <= k and k <= GA.dim and (len(Xk_lis) == GA.bin[k] or len(Xk_lis) == 0), " grade or length wrong ";  # , k, GA.bin[k], len(Xk_lis)  o/p ?? 
    Z = MultiVector(GA.layout); 
    if Xk_lis <> [] : 
      h0 = GA.binsum[k]; 
      for h in range(0, GA.bin[k]) : Z.value[h0+h] = Xk_lis[h]; # end for end if 
    return Z;  # end def 
  
  # Magnitude  ||X||  of multor  X 
  def mag2 (GA, X) : 
    return X.mag2(); # end def 
  
  # Reversion (main anti-automorphism)  X+  of multor  X 
  def rev (GA, X, lo = None, hi = None) :  # local Z; 
    # local i,k,sig,Z; global dim,bin; 
    if lo is None : lo = 0; hi = GA.dim;  # grade options 
    elif hi is None : hi = lo; # end if 
    Z = X.adjoint(); 
    if lo <> 0 or hi <> GA.dim : Z = GA.grades(Z, lo, hi); # end if 
    return Z; # end def 
  
  # Parity transform (main automorphism)  X*  of multor  X 
  def par (GA, X, lo = None, hi = None) : 
    # local i,k,sig,Z; global dim,bin; 
    if lo is None : lo = 0; hi = GA.dim;  # grade options 
    elif hi is None : hi = lo; # end if 
    Z = X.gradeInvol(); 
    if lo <> 0 or hi <> GA.dim : Z = GA.grades(Z, lo, hi); # end if 
    return Z; # end def 
  
  # Dual  Z  of multor  X  (combinatorial: null signatures signed!); 
  def dual (GA, X, lo = None, hi = None) :  # local Z; 
    if lo is None : lo = 0; hi = GA.dim;  # grade options 
    elif hi is None : hi = lo; # end if 
    Z = X.dual(); 
    if lo <> 0 or hi <> GA.dim : Z = GA.grades(Z, lo, hi); # end if 
    return Z;  # end def 
  
  # Add multors: return  Z = X + Y 
  def add (GA, X, Y, lo = None, hi = None) :  # local Z; 
    if lo is None : lo = 0; hi = GA.dim;  # grade options 
    elif hi is None : hi = lo; # end if 
    Z = X.__add__(Y); 
    if lo <> 0 or hi <> GA.dim : Z = GA.grades(Z, lo, hi); # end if 
    return Z; # end def 
  
  # Subtract multors: return  Z = X - Y 
  def sub (GA, X, Y, lo = None, hi = None) :  # local Z; 
    if lo is None : lo = 0; hi = GA.dim;  # grade options 
    elif hi is None : hi = lo; # end if 
    Z = X.__sub__(Y); 
    if lo <> 0 or hi <> GA.dim : Z = GA.grades(Z, lo, hi); # end if 
    return Z; # end def 
  
  # Multiply multors: return  Z = X Y , grades  lo  to  hi  inclusive 
  def mul (GA, X, Y, lo = None, hi = None) :  # local Z; 
    if lo is None : lo = 0; hi = GA.dim;  # grade options 
    elif hi is None : hi = lo; # end if 
    # Z = X.__and__(Y); 
    Z = X.__mul__(Y); 
    if lo <> 0 or hi <> GA.dim : Z = GA.grades(Z, lo, hi); # end if 
    return Z; # end def 
  
  # Dual multiply: return  Z = (X~ Y~)~ , grades  lo  to  hi  inclusive; RK omitted! 
  def lum (GA, X, Y, lo = None, hi = None) : 
    if lo is None : lo = 0; hi = GA.dim;  # grade options 
    elif hi is None : hi = lo; # end if 
    return GA.dual(GA.mul(GA.dual(X), GA.dual(Y), GA.dim-hi, GA.dim-lo));  # end def 
  
  # Left contract multors: return  Z = X _| Y , grades  lo  to  hi  inclusive 
  #   \sum_k,l < <X>_k <Y>_l >_l-k  
  def lecon (GA, X, Y, lo = None, hi = None) :  # local Z; 
    if lo is None : lo = 0; hi = GA.dim;  # grade options 
    elif hi is None : hi = lo; # end if 
    Z = X.lc(Y); 
    if lo <> 0 or hi <> GA.dim : Z = GA.grades(Z, lo, hi); # end if 
    return Z; # end def 
  
  # Right contract multors: return  Z = X _| Y , grades  lo  to  hi  inclusive 
  #   \sum_k,l < <X>_k <Y>_l >_k-l  so  l <= dim-m  &  k = l+m 
  def ricon (GA, X, Y, lo = None, hi = None) : 
    if lo is None : lo = 0; hi = GA.dim;  # grade options 
    elif hi is None : hi = lo; # end if 
    Z = Y.adjoint().lc(X.adjoint()).adjoint();  # patched 
    if lo <> 0 or hi <> GA.dim : Z = GA.grades(Z, lo, hi); # end if 
    return Z; # end def 
  
  # Outer product of multors: return  Z = X ^ Y , grades  lo  to  hi  inclusive 
  #   \sum_k,l < <X>_k <Y>_l >_k+l 
  def wedge (GA, X, Y, lo = None, hi = None) : # local Z; 
    if lo is None : lo = 0; hi = GA.dim;  # grade options 
    elif hi is None : hi = lo; # end if 
    Z = X.__xor__(Y); 
    if lo <> 0 or hi <> GA.dim : Z = GA.grades(Z, lo, hi); # end if 
    return Z; # end def 
  
  # Dual outer product: return  Z = X v Y = (X~ ^ Y~)~ ; RK omitted! 
  def vee (GA, X, Y, lo = None, hi = None) : 
    if lo is None : lo = 0; hi = GA.dim;  # grade options 
    elif hi is None : hi = lo; # end if 
    return GA.dual(GA.wedge(GA.dual(X), GA.dual(Y), GA.dim-hi, GA.dim-lo));  # end def 
  
  # Dorst inner product multors: return  Z = X o Y , grades  lo  to  hi  inclusive 
  #   \sum_k,l < <X>_k <Y>_l >_|k-l| 
  def fat (GA, X, Y, lo = None, hi = None) :  # local Z; 
    if lo is None : lo = 0; hi = GA.dim;  # grade options 
    elif hi is None : hi = lo; # end if 
    Z = (X|Y) + X(0)*Y + X*Y(0) - X(0)*Y(0); 
    if lo <> 0 or hi <> GA.dim : Z = GA.grades(Z, lo, hi); # end if 
    return Z; # end def 
  
  # Hestenes inner product multors: return  Z = X . Y , grades  lo  to  hi  inclusive 
  #   \sum_k,l < <X>_k <Y>_l >_|k-l| 
  def dot (GA, X, Y, lo = None, hi = None) :  # local Z; 
    if lo is None : lo = 0; hi = GA.dim;  # grade options 
    elif hi is None : hi = lo; # end if 
    Z = X.__or__(Y); 
    if lo <> 0 or hi <> GA.dim : Z = GA.grades(Z, lo, hi); # end if 
    return Z; # end def 
  
  # Scalar product multors: return  Z = X * Y 
  #   \sum_k,l < <X>_k <Y>_l >_0 , so  m = 0  &  k = l ; RK omitted! 
  def scalar (GA, X, Y) : # local Z; 
    # local k,l,Z; 
    Z = GA.bld(); 
    for l in range(0, GA.dim+1) : 
      for k in range(0, GA.dim+1) : 
        Z = GA.add(Z, GA.mul(GA.gra(X, k), GA.gra(Y, l), 0)); 
    # end for end for 
    return Z; # end def 
  
  # Conjugation transform of multor: return  Z = (Y+) X Y , grades  lo  to  hi ; RK omitted! 
  def form (GA, X, Y, lo = None, hi = None) : 
    if lo is None : lo = 0; hi = GA.dim;  # grade options 
    elif hi is None : hi = lo; # end if 
    return GA.mul(GA.mul(GA.rev(Y), X), Y, lo, hi);  # end def 
  
  # Sum list of multors ; RK omitted! 
  def addlis (GA, Ylis) :  # local Z,i; 
    Z = GA.bld(); 
    for i in range(0, len(Ylis)) : Z = GA.add(Z, Ylis[i]); # end for 
    return Z; # end def 
  
  # Product list of multors # ; RK omitted! 
  def mullis (GA, Ylis) :  # local Z,i; 
    Z = GA.bld([1]); 
    for i in range(0, len(Ylis)) : Z = GA.mul(Z, Ylis[i]); # end for 
    return Z; # end def 
  
  # Retrieve component of monomial in  X 
  def getcom (GA, X, mon = 0) :  # local k, Xk; 
    assert False, " UNIMPLEMENTED ";  #  --- TEMP ?? 
    assert 0 < mon or mon > GA.len, " illegal monomial "; 
    k = GA.gramon[mon]; 
    return X[k][GA.offmon[mon]] if k >= 0 and X[k] <> [] else 0; # end def 
  
  # Set component of monomial in  X 
  def putcom (GA, X, mon = 0, val = 0) :  # local k; 
    assert False, " UNIMPLEMENTED ";  #  --- TEMP ?? 
    assert 0 < mon or mon > GA.len, " illegal monomial "; 
    k = GA.gramon[mon];  # now build if empty 
    if X[k] == [] : X[k] = [ 0 for i in range(0, GA.bin[k]) ]; # end if 
    X[k][GA.offmon[mon]] = val; 
    return None; # end def 
  
  # Grade  k  and offset  t  of monomial 
  def graset (GA, mon) : 
    assert False, " UNIMPLEMENTED ";  #  --- TEMP ?? 
    assert 0 < mon or mon > GA.len, " illegal monomial "; 
    return GA.gramon[mon], GA.offmon[mon]; # end def 
  
  # Monomial at grade  k  and offset  t 
  def monom (GA, k, t) : 
    assert False, " UNIMPLEMENTED ";  #  --- TEMP ?? 
    return GA.mongo[k][t] if 0 <= k and k <= GA.dim and 0 <= t and t < GA.bin[k]else 0; # end def 
  
  # Distance from multor  X  to zero 			 --- UNUSED ?? 
  def zero (GA, X) : 
    # local Xk, t, err; global GA.eps, GA.len; 
    err = 0.0;  # error 
    Xk = X.__getslice__(0, GA.len); 
    for t in range(0, GA.len) : 
      err = err + abs(Xk[t]); # end for 
    err = err/GA.len; 
    return err; # end def 
  
  # Distance from multors  X, Y  to proportionality 	 --- UNUSED ?? 
  def prop (GA, X, Y) : 
    # local Xk, Yk, Xcs, Ycs, t, err; global GA.eps, GA.len; 
    Xk = X.__getslice__(0, GA.len);  Yk = Y.__getslice__(0, GA.len); 
    Xcs = 0.0;  Ycs = 0.0;  # component sums 
    for t in range(0, GA.bin[k]) : 
      Xcs = Xcs + abs(X[k][t]);  Ycs = Ycs + abs(Y[k][t]); # end for 
    err = 0.0; 
    for t in range(0, GA.bin[k]) : 
      err = err + abs(Xcs*Y[k][t] - Ycs*X[k][t]); # end for 
    err = err/(2.0*Xcs*Ycs) if (Xcs+Ycs)/2.0 > GA.eps else 1.0; 
    return err; # end def 
  
  # Is multor  X  near zero?  Set _eps 				 --- ?? 
  def is_zero (GA, X) : 
    # global _eps; 
    return not(X.__nonzero__()); # end def 
  
  # Are multors  X, Y  nearly proportional?  Set _eps ; fudge! 	 --- ?? 
  def is_prop (GA, X, Y) : 
    # global _eps; 
    return X.__cmp__()(Y); # end def 
  
  # end class 


if False : # if True : # visual check tables 
  GA = multor([+1,+1]); 
  print GA.dim; print GA.len; print GA.pow; print GA.bin; 
  print GA.gensig; print GA.sigsig; print GA.mongo; print GA.gramon; print GA.offmon; 
  print GA.siggo; print GA.doggo; print GA.leng3; print GA.ctrg3; 
  print GA.ilisg3; print GA.jlisg3; print GA.hlisg3; print GA.slisg3; 
# end if 

# Demo integer test of algebraic identities:  Cl(6) , Cl(3,3) ; 
print "clifford wrapper testing:";  # update! 
secs = timeit.default_timer(); 
for (sigs, mags) in \
[ ( [+1,+1,+1,+1,+1,+1], [30517578125, 15625, 12064000] ) \
, ( [+1,+1,+1,-1,-1,-1], [307546875, 3375, -4193280] ) \
#, ( [+1,+1,-1,-1,+0,+0], [530841600, 3600, 38400] ) ] : 
] :  # degenerate unavailable 
  GA = multor(sigs);  # initialise algebra 
  X = GA.mullis([ GA.add(GA.bld([2]), GA.mul(GA.gen(j+1), GA.gen(i+1))) \
  for i in range(0, GA.dim) for j in range(0, i) ]); X; 
  # X = (2 + e1 e2)(2 + e1 e3)(2 + e2 e3) ... 
  Y = GA.mullis([ GA.add(GA.bld([2]), GA.gen(i+1)) 
  for i in range(0, GA.dim) ]); Y;  # Y = (2 + e1)(2 + e2)(2 + e3) ... 
  Z = GA.mullis([ GA.add(GA.bld([2]), GA.mul(GA.bld([i+1]), GA.gen(i+1))) \
  for i in range(0, GA.dim) ]); Z;  # Z = (2 + e1)(2 + 2 e2)(2 + 3 e3) ... 
  
  #assert GA.add(GA.even(Z), GA.odd(Z)) == Z, \
  "clifford:  even(Z) + odd(Z) = Z "; 
  #assert GA.add(GA.grades(Z, 0, 3), GA.grades(Z, 4, 6)) == Z, \
  "clifford:  <Z>_{0,3} + <Z>_{4,6} = Z "; 
  assert [GA.mag2(X), GA.mag2(Y), GA.mag2(Z)] == mags, \
  "clifford:  ||X||, ||Y||, ||Z|| "; 
  assert GA.mag2(X) == GA.lis(GA.mul(GA.rev(X), X))[0], \
  "clifford:  ||X|| = < X+ X >_0 "; 
  assert GA.mul(GA.mul(X, Y), Z) == GA.mul(X, GA.mul(Y, Z)), \
  "clifford:  X (Y Z) = (X Y) Z "; 
  assert GA.mul(X, GA.add(Y, Z)) == GA.add(GA.mul(X, Y), GA.mul(X, Z)), \
  "clifford:  X (Y + Z) = X Y + X Z "; 
  assert GA.lecon(GA.wedge(X, Y), Z) == GA.lecon(X, GA.lecon(Y, Z)), \
  "clifford:  (X ^ Y) _| Z =  X _| (Y _| Z) "; 
  assert GA.add(GA.fat(X, Z), GA.scalar(X, Z)) == GA.add(GA.lecon(X, Z), GA.ricon(X, Z)), \
  "clifford:  (X o Z) + (X * Z) = (X _| Z) + (X |_ Z) "; 
  assert GA.add(GA.fat(X, Z), GA.mul(GA.gra(X), GA.gra(Z))) \
  == GA.add(GA.dot(X, Z), GA.add(GA.mul(GA.gra(X), Z), GA.mul(X, GA.gra(Z)))), \
  "clifford:  (X o Z) + <X>_0 <Z>_0 = (X . Z) + <X>_0 Z  + X <Z>_0 "; 
  assert GA.rev(GA.wedge(Z, Y)) == GA.wedge(GA.rev(Y), GA.rev(Z)), \
  "clifford:  (Z ^ Y)+ = (Y+) ^ (Z+) "; 
  assert GA.wedge(X, GA.wedge(Y, Z)) == GA.wedge(GA.wedge(X, Y), Z), \
  "clifford:  X ^ (Y ^ Z) = (X ^ Y) ^ Z "; 
  assert GA.vee(X, GA.vee(Y, Z)) == GA.vee(GA.vee(X, Y), Z), \
  "clifford:  X v (Y v Z) = (X v Y) v Z "; 
  assert GA.rev(GA.mul(Z, Y)) == GA.mul(GA.rev(Y), GA.rev(Z)), \
  "clifford:  (Z Y)+ = (Y+) (Z+) "; 
  assert GA.par(GA.mul(Z, Y)) == GA.mul(GA.par(Z), GA.par(Y)), \
  "clifford:  (Z Y)* = (Z*) (Y*) "; 
  assert GA.rev(GA.rev(Z)) == Z, \
  "clifford:  (Z+)+ = Z "; 
  assert GA.par(GA.par(Z)) == Z, \
  "clifford:  (Z*)* = Z "; 
  s = 1 - ( GA.dim - sum(GA.sigsig)//2 )%2*2;  # sign change via squared dual 
  assert GA.dual(GA.dual(Z)) == GA.mul(GA.bld([s]), Z), \
  "clifford:  (Z~)~ = +/- Z "; 

# end for 
secs = timeit.default_timer() - secs; 
print " elapsed time in secs ", secs;  # 0.5 sec 

# TODO --- 
# Operators compatible with GAlgebra:  "&", "*" -> "*", "|"  ! 
# Note call func inherited from  __init__  via  X.func() : clumsy! 
# Note RK compos converted to fl. pt. only; note non-degenerate  only!  
# is_zero(), is_prop(): fudged ... 			 --- INCOMPLE ?? 
# monomials:  __getitem__(self, key) , __setitem__(self, key, value) ;  --- INCOMPLE ?? 
# lis(), bld():  __getslice__(self, i, j) , __setslice__(self, i, j, sequ) ?? 
# See clifford-master/test_clifford.py (in /Applications/ ??); 
#   https://media.readthedocs.org/pdf/clifford/latest/clifford.pdf 
#   RK "blade" = grator WFL in pdf doc above! 
# WTF is "projection" --- see L. Dorst [inner.pdf] ?? 
# Tested on GA.givens.py : Cl(4) & Spin_disc() only !! 
# Tested on GA.scripts.py : Apollo, ver. con.  ran OK ; rest N/A 
# User assigning GA.err leaves _eps unchanged: get/put() access ??  --- BUG ?? 
#   from clifford import _eps; print "eps = ", _eps;  #  --- TEMP ?? 
#   in user program outputs  1e-15  always: cannot reset ?! 	 --- FAILS ?? 
#   is_zero(), is_prop() : unusable in practice, substitute WFL versions! 
# GA.even() , GA.odd()  not implemented ?? 

################################################################################

