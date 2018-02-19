
################################################################################

#   "ClifFred" demonstrations for real degenerate Clifford algebras  Cl(p,q,r) . 
# Version 1.2; date: 19/02/18; author: Fred Lunnon <Fred.Lunnon@gmail.com> 
# In command window execute: 
#   cd /Users/fred/fred/python; python -i GA_scripts.py 

# TODOS: optional random fl. pt. inputs? iterative  1/X  demo? 
# Contents list ?? Set timers after initialising? 
# Warm up PRNG using timer? 
# Demo octonions and triality --- to GA_scripts.py ? see  scratch.py 

GAS_version = 1.2;  # update !! 
from GA_multor import *;  # ClifFred 
# from RK_multor import *;  # clifford  --- TEMP ?? 

verbose = True; demons = True;  # print & execute switches, True / False

################################################################################

# Toolkit : random multor generation, etc. 
import timeit; 
from sympy import *;  # symbols, matrices, functions 
import random; 

# Expand symbolic multor X 
def simpex (X) : 
  return [ [ expand(Xki) for Xki in Xk ] for Xk in X ]; # end def  

# Rescale multor to magnitude +/-1 (assumed nonzero) 
def normalise (X) : 
  return GA.mul(X, GA.bld([1.0/sqrt(abs(GA.mag2(X)))])) # end def 

# Random real: component range (-m..+m) ; precision ~ 10E-16 
def rand_unit (m = 1) : 
  return random.uniform(-1, +1); # end def 

# Random vector or grator: grade  k , component range  m  
def rand_grator (k = 1, m = 1) :  # local j; 
  return GA.bld([ rand_unit() for j in range(0, GA.bin[k]) ], k); # end def 

# Random multor: component range  m  
def rand_multor (m = 1) :  # local i; 
  return GA.addlis([ rand_grator(i) for i in range(0, GA.dim+1) ]); # end def 

# Random versor: orthonormal product of  l  vectors, avoiding isotropic 
def rand_versor (l = 2) :  # local i; 
  if l%2 == 0 : X = GA.bld([1]);  
  else : X = rand_grator(); # end if 
  for i in range(0, l//2) : 
    X = GA.mul(X, GA.add(GA.bld([1]), GA.mul(rand_grator(), rand_grator()))); 
    # end for 
  return normalise(X); # end def 

# Commutator product  Z = [X, Y] = X Y - Y X  for multors 
def compro (X, Y) : 
  return GA.sub(GA.mul(X, Y), GA.mul(Y, X));  # end def 

# Jacobi relation:  J = [X, [Y, Z]] + [Y, [Z, X]] + [Z, [X, Y]] 
def jacobi (X, Y, Z) :  # local U, V, W; 
  U = compro(X, compro(Y, Z)); 
  V = compro(Y, compro(Z, X)); 
  W = compro(Z, compro(X, Y)); 
  O = GA.addlis([ U, V, W ]); 
  return O;  # end def 

################################################################################

# Jacobi identity random check for commutator product of multors 
if demons : 
  secs = timeit.default_timer(); 
  print; print "Python/ClifFred/GA_scripts: Jacobi identity, version", GAS_version; print; 
  
  sigs = [+1,+1,+1,+1,+1,+1];  # Cl(6,0,0) 
  #sigs = [+1,+1,+1,-1,-1,-1];  # Cl(3,3,0) 
  #sigs = [+1,+1,-1,-1,0,0];  # Cl(2,2,2) 
  GA = multor(sigs); print "signature ", sigs; print; 
  X = rand_multor();  Y = rand_multor();  Z = rand_multor(); 
  U = jacobi(X, Y, Z);  print "multor ", GA.is_zero(U);  # True 
  X = rand_versor(GA.dim);  Y = rand_versor(GA.dim);  Z = rand_versor(GA.dim); 
  U = jacobi(X, Y, Z);  print "versor ", GA.is_zero(U);  # True 
  
  secs = timeit.default_timer() - secs; 
  print "Elapsed time in secs ", secs;  # 0.47 sec 
# end if 

################################################################################

# Versor conservation random check for outer, contract, inner? product in  Cl(n) 
#   --- must test conservation of generator unit vectors ?? 
#   True True True True True True False  for  Cl(12) x2, Cl(6,6), Cl(4,4,4)  in ~4 Ksec each; 

if demons : 
  secs = timeit.default_timer(); 
  print; print "Python/ClifFred/GA_scripts: versor conserve random, version", GAS_version; print; 
  
  #n = 6; sigs = [+1 for j in range(0, n)];  # Cl(6,0,0) 
  n = 6; sigs = [+1 for j in range(0, (n+1)//2)] \
  + [-1 for j in range(0, (n+0)//2)];  # Cl(3,3,0) 
  #n = 6; sigs = [+1 for j in range(0, (n+2)//3)] \
  #+ [-1 for j in range(0, (n+1)//3)] \
  #+ [+0 for j in range(0, (n+0)//3)];  # Cl(2,2,2) 
  
  GA = multor(sigs); print "signature ", sigs; print; 
  for (k, l) in [ (n-1, n-1), (n-1, n), (n, n-1), (n, n) ] : 
    print "factor grades ", k, l; 
    X = rand_versor(k);  Y = rand_versor(l);  R = rand_versor(1); 
    #T = GA.wedge(X, Y);  W = GA.dot(X, Y);  # (Hestenes) 
    T = GA.wedge(X, Y);  W = GA.fat(X, Y);  # (Dorst) 
    U = GA.ricon(X, Y);  V = GA.lecon(X, Y);  
    print "outer ", GA.is_zero(GA.mul(GA.rev(T), T, 1, GA.dim)), \
    GA.is_zero(GA.form(R, T, 2, GA.dim)); 
    print "right contract ", GA.is_zero(GA.mul(GA.rev(U), U, 1, GA.dim)), \
    GA.is_zero(GA.form(R, U, 2, GA.dim)); 
    print "left contract ", GA.is_zero(GA.mul(GA.rev(V), V, 1, GA.dim)), \
    GA.is_zero(GA.form(R, V, 2, GA.dim)); 
    print "inner fat-dot ", GA.is_zero(GA.mul(GA.rev(W), W, 1, GA.dim)); 
    print; 
  # end for 
  
  secs = timeit.default_timer() - secs; 
  print "Elapsed time in secs ", secs;  # 1.25 sec 
# end if 

################################################################################

# Versor conservation proof for outer, contract, inner? products in  Cl(n) 
#   Grades  (n,n-1)  superfluous since  Y ^ X  =  ((X+) ^ (Y+))+ ; 
#   ricon superfluous since  X |_ Y  =  ( (Y+) _| (X+) )+ ; 

if demons : 
  secs = timeit.default_timer(); 
  print; print "Python/ClifFred/GA_scripts: versor conserve proof, version", GAS_version; print; 
  
  #sigs = [+1,+1,+1,+1];  # spherical  n = 4  (many hrs?) 
  sigs = [+1,+1,+1];  # spherical n = 3  (200 sec) 
  GA = multor(sigs); print "signature ", sigs; print; 
  var("a1, a2, a3, a4, b1, b2, b3, b4, c1, c2, c3, c4, d1, d2, d3, d4, \
  e1, e2, e3, e4, f1, f2, f3, f4, g1, g2, g3, g4, h1, h2, h3, h4"); 
  A = GA.bld([a1, a2, a3], 1); # , a4 
  B = GA.bld([b1, b2, b3], 1); # , b4 
  C = GA.bld([c1, c2, c3], 1); # , c4 
  D = GA.bld([d1, d2, d3], 1); # , d4 
  E = GA.bld([e1, e2, e3], 1); # , e4 
  F = GA.bld([f1, f2, f3], 1); # , f4 
  G = GA.bld([g1, g2, g3], 1); # , g4 
  H = GA.bld([h1, h2, h3], 1); # , h4 
  Xfac = [A, B, C, D]; Yfac = [E, F, G, H]; 
  
  #for (k, l) in [ (2, 2), (2, 3), (3, 2), (3, 3), (3, 4), (4, 3), (4, 4)] : 
  for (k, l) in [ (2, 2), (2, 3), (3, 2), (3, 3)] : 
    print "factor grades ", k, l; 
    X = GA.mullis([Xfac[i] for i in range(0, k)]); 
    Y = GA.mullis([Yfac[i] for i in range(0, l)]); 
    U = GA.wedge(X, Y);  V = GA.lecon(X, Y);  W = GA.fat(X, Y);  
    print simpex(GA.mul(GA.rev(U), U, 1, GA.dim)); 
    print simpex(GA.mul(GA.rev(V), V, 1, GA.dim)); 
    print simpex(GA.mul(GA.rev(W), W, 1, GA.dim)); 
    print;  # scalar, scalar, not scalar when  n = 4 ! 
  # end for 
  
  secs = timeit.default_timer() - secs; 
  print "Elapsed time in secs ", secs;  # 195 sec 
# end if 

################################################################################

# Symbolic proof of inverse formula for rank  n = 4  --- 
# if  Z := X X+  then  <Z>_2 = <Z>_3 = 0 ; 
# if  W := Z (2 <Z>_0 - Z)  then  <W>_1 = <W>_2 = <W>_3 = <W>_4 = 0 ; 
# hence  X Y = c  scalar where  Y = X+ (2 <X X+>_0 - X X+) . 

# Expand symbolic multor X 
def simpex (X) : 
  return [ [ expand(Xki) for Xki in Xk ] for Xk in X ]; # end def  

if demons : 
  secs = timeit.default_timer(); 
  print; print "Python/ClifFred/GA_scripts: inverse formula proof, version", GAS_version; print; 
  
  GA = multor([+1,+1,+1,+1]);  #  n = 4  dummy 
  var("x, x1, x2, x3, x4, x12, x13, x23, x14, x24, x34, x123, x124, x134, x234, x1234"); 
  X = GA.addlis([ GA.bld([x], 0), GA.bld([x1, x2, x3, x4], 1), \
  GA.bld([x12, x13, x23, x14, x24, x34], 2), GA.bld([x123, x124, x134, x234], 3), \
  GA.bld([x1234], 4)]); 
  Z = GA.addlis([ GA.bld([x], 0), GA.bld([x1, x2, x3, x4], 1), GA.bld([x1234], 4)]); 
  
  for siglis in [ \
  [+1,+1,+1,+1], [+1,+1,+1,-1], [+1,+1,-1,-1], [+1,-1,-1,-1], [-1,-1,-1,-1], \
  [+1,+1,+1,0], [+1,+1,-1,0], [+1,-1,-1,0], [-1,-1,-1,0], [+1,+1,0,0], \
  [+1,-1,0,0], [-1,-1,0,0], [+1,0,0,0], [-1,0,0,0], [0,0,0,0] ] :  # 15 cases 
    print "signature ", siglis; 
    GA = multor(siglis);  # initialise algebra 
    print "X X+ = ", simpex(GA.mul(X, GA.rev(X)));  # grades 0,1,4  
    print "Z (2 <Z>_0 - Z) = ", simpex(GA.mul(Z, GA.sub(GA.bld([2*x]), Z)));  # grade 0 
  # end for 
  
  secs = timeit.default_timer() - secs; 
  print "Elapsed time in secs ", secs;  # 0.53 sec 
# end if 

################################################################################

# Octonions and Moufang identities: 
#   triality remains unimplemented since Lounesto sect. 9 apparently incorrect! 

# Pertti Lounesto "Octonions and Triality" 
# Advances in Applied Clifford Algebras vol. 11 no.2, 191--213 (2001) , 
# sect. 3 : octonion product for paravectors in  Cl(0,7)  via 
#   A o B  ==  < A B(1 + W)(1 - e_1234567) >_0..1 ; 
#   A o B  ==  < A B(1 - V) >_{0,1} ; 
# for multors in  Cl(0,8)  via 
#   A o B  ==  < A e_8 B(1 + W)(1 - J) >_1  for multors  A,B , identity == e_8 ; 

# Convert  Cl(0,7)  paravector  X  to/from  Cl(8)^0  multor  Y  in  Cl(0,8) : 
#   assumes  X  omits  e_8 ,  Y  is even; 
def con7to8(X) :  # local me8; 
  me8 = GA.sub(GA.bld(), GA.gen(8));  # - e_8 
  return GA.sub(GA.even(X), GA.mul(GA.odd(X), me8));  # end def 

def con8to7(Y) :  # local me8, Xe, Xo; 
  me8 = GA.sub(GA.bld(), GA.gen(8));  # - e_8 
  Xe = GA.mul(GA.wedge(Y, me8), me8);  # - even(Y) 
  Xo = GA.mul(GA.add(Y, Xe), me8);  # odd(Y) 
  return GA.sub(Xo, Xe);  # end def 

# Octonion product of paravectors in  Cl(0,n)  for  n >= 7 : 
#   < X Y (1 - V) >_{0,1} , with  V = e_124 + ... + e_713 ; 
#   note non-associative & non-commutative 
def omp7(X, Y) :  # local V, Z, fano, vert, line; 
  fano = [ [1,2,4], [2,3,5], [3,4,6], [4,5,7], [5,6,1], [6,7,2], [7,1,3] ]; 
  V = GA.addlis([ GA.mullis( [GA.gen(vert) for vert in line]) for line in fano ]); 
  Z = GA.mul(GA.mul(X, Y), GA.sub(GA.bld([1]), V)); 
  return GA.add(GA.gra(Z, 0), GA.gra(Z, 1));  # end def 

# Octonion product of multors in  Cl(0,n)  for  n >= 8 : safe version 
def omp8(X, Y) : 
  return con7to8(omp7(con8to7(X), con8to7(Y)));  # end def 


# Prove Moufang identities for symbolic octonions in  Cl(0,7) 
#   ((X Y)X)Z = X(Y(X Z);  Z((X Y)X) = ((Z X)Y)X;  (X(Y Z))X = (X Y)(Z X);  
verbose = True; demons = True;  # print & execute switches, True / False
if demons : 
  secs = timeit.default_timer(); 
  print; print "Python/ClifFred/GA_scripts: octonions, Moufang, triality :  ", GAS_version; print; 
  
  var("x0, x1, x2, x3, x4, x5, x6, x7, x8, y0, y1, y2, y3, y4, y5, y6, y7, y8, z0, z1, z2, z3, z4, z5, z6, z7, z8"); 
  n = 7; sigs = [-1 for j in range(0, n)];  #  Cl(0,7) 
  GA = ClifFred(sigs); print "signature ", sigs; print; 
  
  X = GA.add(GA.bld([x0], 0), GA.bld([x1, x2, x3, x4, x5, x6, x7], 1)); 
  Y = GA.add(GA.bld([y0], 0), GA.bld([y1, y2, y3, y4, y5, y6, y7], 1)); 
  Z = GA.add(GA.bld([z0], 0), GA.bld([z1, z2, z3, z4, z5, z6, z7], 1)); 
  simpex(GA.sub( omp7(omp7(X, omp7(Y, Z)), X), omp7(omp7(X, Y), omp7(Z, X)) ));  
  # (X(Y Z))X = (X Y)(Z X) 
  simpex(GA.sub( omp7(omp7(omp7(X, Y), X), Z), omp7(X, omp7(Y, omp7(X, Z))) )); 
   # ((X Y)X)Z = X(Y(X Z) 
  simpex(GA.sub( omp7(Z, omp7(omp7(X, Y), X)), omp7(omp7(omp7(Z, X), Y), X) )); 
   # Z((X Y)X) = ((Z X)Y)X 
  
  secs = timeit.default_timer() - secs; 
  print "Elapsed time in secs ", secs;  # 16.5 sec 
# end if 

################################################################################

# Regular tetrahedron in Euclidean / projective geometry --- 
#   generators  x_,y_,z_,o_  represent coordinate & infinity planes; 
#   vector  a x_ + b y_ + c z_ + d o_ , plane  a x + b y + c z = -d ; 
#   covector  s xyz_ + r xyo_ - q xzo_ + p yzo_ , point  (p/s, q/s, r/s) ; 
# Todos: rescale versors by 1/GCD; 
#   investigate false bi-alt.  K = GA.lum(L, M, 2); 

if demons : 
  secs = timeit.default_timer(); 
  print; print "Python/ClifFred/GA_scripts: Tetrahedron, version", GAS_version; print; 
  
  GA = multor([+1,+1,+1,+0]);  # 3-space Euclidean geometry  Cl(3,0,1) 
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
  
  secs = timeit.default_timer() - secs; 
  print "Elapsed time in secs ", secs;  # 0.0057 sec 
# end if 

################################################################################

# Apollonian problem in Lie-sphere 2-space: output both cycles 
#   (directed circles) anti-tangent to three input cycles.  

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
  X1 = GA.lis(X, 1);  # extract vector 
  xx = [ X1[i] for i in range(0, n) ]; r = X1[n+2]; 
  o = (X1[n] - X1[n+1])/2.0; p = (X1[n] + X1[n+1])/2.0;  # fl. pt. divide! 
  if abs(p) > CPR*abs(o) and abs(p) > CPR*abs(r) : 
    if abs(r) > CPR*abs(2*p) : return [ 1 ] + [ X1[i]/(2.0*p) for i in range(0, n) ] + [ r/(2.0*p) ]; # sphere 
    else : return [ 1 ] + [ X1[i]/(2.0*p) for i in range(0, n) ] + [ r/(2.0*p) ]; # end if : point 
  elif abs(r) > CPR*abs(o) : return [ o/(1.0*r) ] + [ X1[i]/(1.0*r) for i in range(0, n) ] + [ oo ];  # prime 
  elif o <> 0 : return [ 1 ] + [ 0 for i in range(0, n) ] + [ oo ];  # infinity :  abs(o) > CPR  ?? 
  else : return [ 0 ] + [ 0 for i in range(0, n) ] + [ 0 ]; # end if : zero 
  # end def 

# Lie product (X | Y) = -1/2(tang dist)^2, cos angle - 1, if normalised. 
# Generalise to same-parity versors via inn(X, Y) = <X+ Y>_0 ? 
def inn (X, Y) :  # --- UNUSED ?? 
  return GA.lis(GA.mul(X, Y, 0))[0]; # end def 

# Lie cross-ratio of 4 cycles --- treat cases oo and 0/0 ?? 
def crorat (X, Y, Z, W) :  # --- UNUSED ?? 
  return (inn(X, Y)*inn(Z, W))/(inn(X, Z)*inn(Y, W)); # end def 

# Vector Z with given base cycle Y, angle s, tangent length t: 
#   Z = [Z_x,...,Z_r] where Z_i = Y_i for i <> r, Z_r = Y_r cos s where s equals 
#   angle made by the fixed cycles X <> Y of Z with Y; 
#   Z_i = Y_i for i <> o, Z_o = Y_o + (t Y_p)^2 where t equals tangent length 
#   to the fixed cycles X <> Y of Z from Y. 
def modvec_t2s2 (X, tansqd, sin2son2) : # t^2, sin^2(s/2); not both nonzero! 
  # local X1,i; 
  X1 = GA.lis(X, 1); # extract vector 
  return GA.bld([ X1[i] for i in range(0, n) ] + [ X1[n]+tansqd/2, X1[n+1]-tansqd/2, X1[n+2]*(1-2*sin2son2) ], 1); # end def 

# Common tangents (assumed real) of  n+1  cycles or vectors in  n-space; 
#   odd permutation of inputs reverses order of outputs. 
# Eigencycle algorithm: assuming ||<Z>_2|| = -c^2 <= 0, then c +/- <Z>_2 singular; 
#   so (c +/- <Z>_2) W (c -/+ <Z>_2) = X,Y gives the eigencycles for any test 
#   cycle W, say W = oo. What happens if rotation Z is parabolic tangential? 
# Maybe use -sqrt(-magX) or +sqrt(-magX) according to "sign" of X [e.g. of sum 
#   of comps], to establish continuity of distinct tangents? 
def com_tang (Ylis) : 
  # local i, X, X1, X2, X11, X21, magX, V, pimin2; 
  
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
if demons : 
  secs = timeit.default_timer(); 
  print; print "Python/ClifFred/GA_scripts: Apollonian, version", GAS_version; print; 
  
  GA = multor([+1,+1,+1,-1,-1]);  # 2-space Lie-sphere geometry  Cl(3,2) 
  n = GA.dim-3; oo = 1.0E99;  # --- TEMP ?? 
  
  # cycle format  Xc = [1, x, y, r]  denotes  X  at centre  (x,y) , radius  r ; 
  Y1c = [1, 15, 0, 15]; Z1 = modvec_t2s2(lie(Y1c), 0, 1); 
  Y2c = [1, -15, 0, 15]; Z2 = modvec_t2s2(lie(Y2c), 0, 1); 
  Y3c = [1, 0, 20, 10]; Z3 = modvec_t2s2(lie(Y3c), 0, 1); 
  if verbose : print "Y1c, Z1, Y2c, Z2, Y3c, Z3 "; print Y1c; print Z1; print Y2c; print Z2; print Y3c; print Z3; print; 
  
  X1, X2 = com_tang([Z1, Z2, Z3]); 
  X1m = sqrt(abs(GA.mag2(X1))); 
  X2m = sqrt(abs(GA.mag2(X2))); # isotropic vectors? 
  X1c = cyc(X1); X2c = cyc(X2); # extract [1, x, y, r] 
  X1c, X2c;  # X1c := [1, 0, 8, 2]; X2c := [1, 0, 0, -30]; 
  if verbose : print "X1c, X1m, X1, X2c, X2m, X2 "; print X1c, X1m; print X1; print X2c, X2m; print X2; print; 
  
  secs = timeit.default_timer() - secs; 
  print "Elapsed time in secs ", secs;  # 0.018 sec 
# end if 

################################################################################

