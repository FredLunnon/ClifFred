
################################################################################

#   "ClifFred" demonstrations for real degenerate Clifford algebras  Cl(p,q,r) . 
# Version 1.1; date: 20/08/16; author: Fred Lunnon <Fred.Lunnon@gmail.com> 
# In command window execute: 
#   cd /Users/fred/fred/python; python -i GA_scripts.py 
# Todos: optional random fl. pt. inputs? 

GAS_version = 1.1;  # update !! 

################################################################################

# Regular tetrahedron in Euclidean / projective geometry --- 
#   generators  x_,y_,z_,o_  represent coordinate & infinity planes; 
#   vector  a x_ + b y_ + c z_ + d o_ , plane  a x + b y + c z = -d ; 
#   covector  s xyz_ + r xyo_ - q xzo_ + p yzo_ , point  (p/s, q/s, r/s) ; 
# Todos: rescale versors by 1/GCD; 
#   investigate false bi-alt.  K = GA.lum(L, M, 2); 
import timeit; 
from GA_multor import *;  # ClifFred 

verbose = True; demons = True;  # print & execute switches, True / False
secs = timeit.default_timer(); 
if demons : 
  print; print "Python/ClifFred/GA_scripts: Tetrahedron, version", GAS_version; print; 
  
  GA = ClifFred([+1,+1,+1,+0]);  # 3-space Euclidean geometry  Cl(3,0,1) 
  x_ = GA.gen(1); y_ = GA.gen(2); z_ = GA.gen(3); o_ = GA.gen(4); 
  print "  coordinate planes and oo = ", x_, y_, z_, o_; print; 
  
  FF = [ GA.bld([+1, +1, +1, 1], 1), \
  GA.bld([+1, -1, -1, 1], 1), \
  GA.bld([-1, +1, -1, 1], 1), \
  GA.bld([-1, -1, +1, 1], 1) ]; 
  print "  faces F_i = ", FF; print;  # 4 face planes 
  
  PP = [ GA.gra(GA.mullis([ FF[i], FF[j], FF[k] ]), 3) \
  for (i,j,k) in [ (1,2,3), (2,3,0), (3,0,1), (0,1,2) ] ]; 
  print "  vertices P_j = ", PP; print;  # vertex point meets 3 adjacent faces 
  
  FdotP = [ GA.lis(GA.lum(F, P, 0, 0))[0] for P in PP for F in FF ]; 
  print "  F_i . P_j = ", FdotP; print;  #  F_i meets P_j  for  i <> j 
  
  FFdup = [ GA.lum(GA.lum(PP[i], PP[j]), PP[k], 1, 1) \
  for (i,j,k) in [ (1,2,3), (2,3,0), (3,0,1), (0,1,2) ] ]; 
  print "  F_i (scaled) = ", FFdup; print;  # face joins 3 adjacent vertices 
  
  # bivector blade  L = [ P, F ] represents (Pluecker) line: 
  # P = <L \o>_3  direction point at infinity, 
  # F = < (L~ \o)~ >_1  moment plane joining line to origin. 
  LL = [ GA.mul(FF[i], FF[j], 2) \
  for (i,j) in [ (0,1), (0,2), (1,2), (0,3), (1,3), (2,3) ] ];  
  print "  edges L_ij = ", LL; print;  # edge line = meet of 2 adjacent faces 
  
  # Grassmann relation: lines are versors 
  LdagL = [ GA.mul(L, GA.rev(L)) for L in LL ]; 
  print "  (L_i+) L_i = ", LdagL; print;  #  (L+) L  is scalar 
  
  # edge-lengths squared 
  eleng2s = [ GA.mag2(GA.lum(P, Q, 2)) /1.0**2 / GA.mag2(GA.mul(P, Q, 0)) \
  for P in PP for Q in PP ]; 
  print "  edge-lengths squared = ", eleng2s; print;  # 8 
  
  # face-areas squared (divide by 4 ?) 
  farea2s = [ GA.mag2(GA.lum(GA.lum(PP[i], PP[j]), PP[k], 1)) \
  /2.0**2 / GA.mag2(GA.mul(GA.mul(PP[i], PP[j]), PP[k], 3)) \
  for (i,j,k) in [ (1,2,3), (2,3,0), (3,0,1), (0,1,2) ] ]; 
  print "  face-areas squared = ", farea2s; print;  # 12 
  
  # volume (omitted sign?) 
  vol = GA.lis(GA.lum(GA.lum(GA.lum(PP[0], PP[1]), PP[2]), PP[3], 0))[0] \
  /6.0 / GA.lis(GA.mul(GA.mul(GA.mul(PP[0], PP[1]), PP[2]), PP[3], 0))[0]; 
  print "  volume = ", vol; print;  # 8/3 
  
  # altitude line  N = < F P >_2 , foot point  Q = < F M >_3 , length squared = ||< (F~ P~)~ >_0|| / ||< F P >_2|| 
  F = FF[0]; P = PP[0]; 
  N = GA.mul(F, P, 2);  Q = GA.mul(F, N, 3); 
  FNc = GA.mag2(GA.mul(F, N, 1));  # 0 , so  N  perp'  F  
  len2 = GA.mag2(GA.lum(F, P, 0)) /1.0**2 / GA.mag2(GA.mul(F, P, 2)); 
  print "  altitude line, foot point, cosine, length squared = ", N, Q, FNc, len2; print;  # 16/3 
  
  # bi-altitude line  N = < L M >_2 , meets  L, M  perp'ly, length squared = ||< (L~ M~)~ >_0|| / ||< L M >_2|| 
  L = LL[0]; M = LL[5];  # opposite edges 
  N = GA.mul(L, M, 2); N;  #  N = line  y = z = 0 
  LNc = GA.mag2(GA.mul(L, N, 0)); LNd = GA.mag2(GA.lum(L, N, 0));  # 0, 0 
  MNc = GA.mag2(GA.mul(M, N, 0)); MNd = GA.mag2(GA.lum(M, N, 0));  # 0, 0 
  len2 = GA.mag2(GA.lum(L, M, 0)) /1.0**2 / GA.mag2(GA.mul(L, M, 2)); 
  print "  bi-altitude line, cosines, dists, length squared = ", N, LNc, MNc, LNd, MNd, len2; print;  # 4 

# end if 
secs = timeit.default_timer() - secs; 
print "Elapsed time in secs ", secs; 

################################################################################

# Apollonian problem in Lie-sphere 2-space: output both cycles 
#   (directed circles) anti-tangent to three input cycles.  
import timeit; 
from GA_multor import *;  # ClifFred 
from sympy import *;  # matrices & functions 

# Lie coordinates X of cycle in n-space (version 7); 
# homog cycle sphere/point [2p, a, b, r], prime [o, a, b, oo] 
def lie (wxyr) : # local xx,i,r,o,p,sumx2; global n,oo; 
  p = wxyr[0]/2.0; xx = [ wxyr[i+1] for i in range(0, n) ]; r = wxyr[n+1]; 
  sumx2 = 0; 
  for i in range(0, n) : sumx2 = sumx2 + xx[i]*xx[i]; 
  o = wxyr[0] if r == oo else -(sumx2 - r*r)/2.0; # end if 
  return GA.bld( xx + [ o, -o, sqrt(sumx2) ] if r == oo else [ 0 for i in range(0, n) ] + [ +1, -1, 0 ] if p == 0 else xx + [ p + o/(2.0*p), p - o/(2.0*p), r ], 1); # end if 
  # prime, infinity, sphere resp. 
  # end def 

# Homog cycle from Lie coordinates X --- assumes ||X|| = 0, normalising; 
#   extreme spheres rounded to points or lines or oo. 
# Mods: adjust rounded prime so distance from origin is unaltered ?? 
#   Round small p,r to infinity; small p,r,o to zero ?? 
#   Test for approx. cycle; round to "nearest" cycle ?? 
CPR = 0.000001; # line & point cutoff power & radius  --- use GA.eps ??  *** 
def cyc (X) : # local X1,p,o,xx,i,r; global CPR,n,oo; 
  X1 = X[1];  # extract vector 
  xx = [ X1[i] for i in range(0, n) ]; r = X1[n+2]; 
  o = (X1[n] - X1[n+1])/2.0; p = (X1[n] + X1[n+1])/2.0;  # fl. pt. divide! 
  if abs(p) > CPR*abs(o) and abs(p) > CPR*abs(r) : 
    if abs(r) > CPR*abs(2*p) : return [ 1 ] + [ X1[i]/(2.0*p) for i in range(0, n) ] + [ r/(2.0*p) ]; # sphere 
    else : return [ 1 ] + [ X1[i]/(2.0*p) for i in range(0, n) ] + [ r/(2.0*p) ]; # end if : point 
  elif abs(r) > CPR*abs(o) : return [ o/(1.0*r) ] + [ X1[i]/(1.0*r) for i in range(0, n) ] + [ oo ];  # prime 
  elif o <> 0 : return [ 1 ] + [ 0 for i in range(0, n) ] + [ oo ];  # infinity :  abs(o) > CPR  ?? 
  else : return [ 0 ] + [ 0 for i in range(0, n) ] + [ 0 ]; # end if : zero 
  # end def 

# Lie product (X | Y) = -1/2(tang dist)^2, cos angle - 1, normalised. 
def inn (X, Y) : 
  return GA.mul(GA.rev(X), Y, 0, 0)[0][0]; # end def 

def inn0 (X, Y) :  # local j; 
  return sum([ X[1][j]*Y[1][j]*GA.gensig[j] for j in range(0, GA.dim) ]); # end def 

# Lie cross-ratio of 4 cycles --- treat cases oo and 0/0 ?? 
def crorat (X, Y, Z, W) : 
  return (inn(X, Y)*inn(Z, W))/(inn(X, Z)*inn(Y, W)); # end def 
  # Generalise to same-parity versors via inn(X, Y) = <X^+ Y>_0 ?? 

# Vector Z with given base cycle Y, angle s, tangent length t: 
#   Z = [Z_x,...,Z_r] where Z_i = Y_i for i <> r, Z_r = Y_r cos s where s equals 
#   angle made by the fixed cycles X <> Y of Z with Y; 
#   Z_i = Y_i for i <> o, Z_o = Y_o + (t Y_p)^2 where t equals tangent length 
#   to the fixed cycles X <> Y of Z from Y. 
def modvec_t2s2 (X, tansqd, sin2son2) : # t^2, sin^2(s/2); not both nonzero! 
  # local X1,i; 
  X1 = X[1]; # extract vector 
  return GA.bld([ X1[i] for i in range(0, n) ] + [ X1[n]+tansqd/2, X1[n+1]-tansqd/2, X1[n+2]*(1-2*sin2son2) ], 1); # end def 

# Common tangents (assumed real) of  n+1  cycles or vectors in  n-space; 
#   odd permutation of inputs reverses order of outputs. 
# Eigencycle algorithm: assuming ||<Z>_2|| = -c^2 <= 0, then c +/- <Z>_2 singular; 
#   so (c +/- <Z>_2) W (c -/+ <Z>_2) = X,Y gives the eigencycles for any test 
#   cycle W, say W = oo. What happens if rotation Z is parabolic tangential? 
# Maybe use -sqrt(-magX) or +sqrt(-magX) according to "sign" of X [e.g. of sum 
#   of comps], to establish continuity of distinct tangents? 
def com_tang (Ylis) : 
  # local i, X, X1, X2, magX, V, pimin2; 
  
  # Extract eigencycles from bivector dual of product 
  X = GA.mullis(Ylis); 
  # X = GA.mul(X, GA.J, 2);  # pseudar dual? 
  X = GA.gra(GA.dual(X), 2);  # reflect dual? 
  magX = min(0, GA.mag2(X)); # check neg, else solutions complex? 
  X = GA.add(GA.bld([sqrt(-magX)], 0), X); # isotropic 
  pimin2 = 1.1415926536;  # "random" cycle 
  V = lie([pimin2**i for i in range(0, GA.dim-1)]); 
  X1 = GA.form(V, X, 1, 1); # isotropic solutions 
  X2 = GA.form(V, GA.rev(X), 1, 1); 
  return X1, X2; # end def 

# Construct antitangent quartets, curvatures = [-1,2,2,3]/30  etc. 
#   --- cycle norms = O(1.0E-5) , yet component error = O(1.0E-14) ! 
verbose = True; demons = True;  # print & execute switches, True / False
secs = timeit.default_timer(); 
if demons : 
  print; print "Python/ClifFred/GA_scripts: Apollonian, version", GAS_version; print; 
  
  GA = ClifFred([+1,+1,+1,-1,-1]);  # 2-space Lie-sphere geometry  Cl(3,2) 
  n = GA.dim-3; oo = 1.0E99;  # --- TEMP ?? 
  
  # cycle format  Xc = [1, x, y, r]  denotes  X  at centre  (x,y) , radius  r ; 
  Y1c = [1, 15, 0, 15]; Z1 = modvec_t2s2(lie(Y1c), 0, 1); 
  Y2c = [1, -15, 0, 15]; Z2 = modvec_t2s2(lie(Y2c), 0, 1); 
  Y3c = [1, 0, 20, 10]; Z3 = modvec_t2s2(lie(Y3c), 0, 1); 
  if verbose : print "Y1c, Z1, Y2c, Z2, Y3c, Z3 "; print Y1c; print Z1; print Y2c; print Z2; print Y3c; print Z3; print; 
  
  X1, X2 = com_tang([Z1, Z2, Z3]); 
  X1m = sqrt(abs(GA.mag2(X1)))/(X1[1][n]+X1[1][n+1]); 
  X2m = sqrt(abs(GA.mag2(X2)))/(X2[1][n]+X2[1][n+1]); 
  X1c = cyc(X1); X2c = cyc(X2); 
  X1c, X2c;  # X1c := [1, 0, 8, 2]; X2c := [1, 0, 0, -30]; 
  if verbose : print "X1c, X1m, X1, X2c, X2m, X2 "; print X1c, X1m; print X1; print X2c, X2m; print X2; print; 
  
# end if 
secs = timeit.default_timer() - secs; 
print "Elapsed time in secs ", secs; 

################################################################################

# Versor conservation random check for outer, contract, inner? product in  Cl(n) 
import random; 
import timeit; 
from GA_multor import *;  # ClifFred 
from sympy import *;  # symbols, matrices, functions 

# Expand symbolic multor X 
def simul (X) : 
  return [ [ expand(Xki) for Xki in Xk ] for Xk in X ]; # end def  

def normalise (X) : 
  return GA.mul(X, GA.bld([1.0/sqrt(abs(GA.mag2(X)))])) # end def 

# Random real: component range (-1..+1) ; precision ~ 10E-16 
def rand_unit () : 
  return random.uniform(-1, +1); # end def 

# Random real list 
def rand_list () :  # local j, ranlis, mag; 
  return [ rand_unit() for j in range(0, GA.dim) ];  # end def 

# Random orthonormal  Cl(p,q,r)  product of  l  vectors, avoiding isotropic 
def rand_versor (l) :  # local i; 
  if l%2 == 0 : X = GA.bld([1]);  
  else : X = GA.bld(rand_list(), 1); # end if 
  for i in range(0, floor(l/2)) : 
    X = GA.mul(X, GA.add(GA.bld([1]), GA.mul(GA.bld(rand_list(), 1), 
      GA.bld(rand_list(), 1)))); # end for 
  return normalise(X); # end def 

# Random test for outer, contract, inner in  Cl(n)  etc. : 
#   --- must test conservation of generator unit vectors ?? 
#   True True True True True True False  for  Cl(12) x2, Cl(6,6), Cl(4,4,4)  in ~4 Ksec each; 

verbose = True; demons = True;  # print & execute switches, True / False
secs = timeit.default_timer(); 
if demons : 
  print; print "Python/ClifFred/GA_scripts: versor conserve random, version", GAS_version; print; 
  
  #n = 12; sigs = [+1 for j in range(0, n)]; 
  #n = 12; sigs = [+1 for j in range(0, (n+1)//2)] + [-1 for j in range(0, (n+0)//2)]; 
  #n = 12; sigs = [+1 for j in range(0, (n+2)//3)] + [-1 for j in range(0, (n+1)//3)] + [+0 for j in range(0, (n+0)//3)]; 
  n = 8; sigs = [+1 for j in range(0, (n+2)//3)] + [-1 for j in range(0, (n+1)//3)] + [+0 for j in range(0, (n+0)//3)]; 

  GA = ClifFred(sigs); print "signature ", sigs; print; 
  for (k, l) in [ (n-1, n-1), (n-1, n), (n, n-1), (n, n) ] : 
    print "factor grades ", k, l; 
    X = rand_versor(k);  Y = rand_versor(l);  R = rand_versor(1); 
    T = GA.wedge(X, Y);  W = GA.dot(X, Y);  # (naive) 
    U = GA.ricon(X, Y);  V = GA.lecon(X, Y);  
    print "outer ", GA.zero(GA.mul(GA.rev(T), T, 1, GA.dim)), \
    GA.zero(GA.form(R, T, 2, GA.dim)); 
    print "right contract ", GA.zero(GA.mul(GA.rev(U), U, 1, GA.dim)), \
    GA.zero(GA.form(R, U, 2, GA.dim)); 
    print "left contract ", GA.zero(GA.mul(GA.rev(V), V, 1, GA.dim)), \
    GA.zero(GA.form(R, V, 2, GA.dim)); 
    print "inner fat-dot ", GA.zero(GA.mul(GA.rev(W), W, 1, GA.dim)); 
    print; 
  # end for 

# end if 
secs = timeit.default_timer() - secs; 
print "Elapsed time in secs ", secs; 

################################################################################

# Versor conservation proof for outer, contract, inner? products in  Cl(n) 
#   Grades  (n,n-1)  superfluous since  Y ^ X  =  ((X+) ^ (Y+))+ ; 
#   ricon superfluous since  X |_ Y  =  ( (Y+) _| (X+) )+ ; 
import timeit; 
from GA_multor import *;  # ClifFred 
from sympy import *;  # symbols, matrices, functions 

# Expand symbolic multor X 
def simul (X) : 
  return [ [ expand(Xki) for Xki in Xk ] for Xk in X ]; # end def  
  
verbose = True; demons = True;  # print & execute switches, True / False
secs = timeit.default_timer(); 
if demons : 
  print; print "Python/ClifFred/GA_scripts: versor conserve proof, version", GAS_version; print; 
  
  #sigs = [+1,+1,+1,+1];  # spherical  n = 4  (many hrs?) 
  sigs = [+1,+1,+1];  # spherical n = 3  (200 sec) 
  GA = ClifFred(sigs); print "signature ", sigs; print; 
  var("a1, a2, a3, a4, b1, b2, b3, b4, c1, c2, c3, c4, d1, d2, d3, d4, \
  e1, e2, e3, e4, f1, f2, f3, f4, g1, g2, g3, g4, h1, h2, h3, h4"); 
  A = GA.bld([a1, a2, a3, a4], 1); 
  B = GA.bld([b1, b2, b3, b4], 1); 
  C = GA.bld([c1, c2, c3, c4], 1); 
  D = GA.bld([d1, d2, d3, d4], 1); 
  E = GA.bld([e1, e2, e3, e4], 1); 
  F = GA.bld([f1, f2, f3, f4], 1); 
  G = GA.bld([g1, g2, g3, g4], 1); 
  H = GA.bld([h1, h2, h3, h4], 1); 
  Xfac = [A, B, C, D]; Yfac = [E, F, G, H]; 
  
  #for (k, l) in [ (2, 2), (2, 3), (3, 2), (3, 3), (3, 4), (4, 3), (4, 4)] : 
  for (k, l) in [ (2, 2), (2, 3), (3, 2), (3, 3)] : 
    print "factor grades ", k, l; 
    X = GA.mullis([Xfac[i] for i in range(0, k)]); 
    Y = GA.mullis([Yfac[i] for i in range(0, l)]); 
    U = GA.wedge(X, Y);  V = GA.lecon(X, Y);  W = GA.dot(X, Y);  
    print simul(GA.mul(GA.rev(U), U, 1, GA.dim)); 
    print simul(GA.mul(GA.rev(V), V, 1, GA.dim)); 
    print simul(GA.mul(GA.rev(W), W, 1, GA.dim)); 
    print;  # scalar, scalar, not scalar when  n = 4 ! 
  # end for 
  
# end if 
secs = timeit.default_timer() - secs; 
print "Elapsed time in secs ", secs; 

################################################################################

# Demo: invert() ; versor conservation ?? 

