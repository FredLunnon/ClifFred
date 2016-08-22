
################################################################################

#   Python/SymPy/GA_multor program source for numerical frame transformation 
# and versor decomposition into fixed-axis (Givens') rotations, with demo & test 
# harness, for real Clifford algebra  Cl(p,q,r) . 
# Version 3.1; date: 30/07/16; author: Fred Lunnon <Fred.Lunnon@gmail.com> 

# In command window execute: 
#     cd /Users/fred/fred/python; python -i frames_givens_demo.py 
# or
#     cd /Users/fred/fred/python; python 
#     from frames_givens_demo import *; 

from sympy import *;  # matrices & functions; for GAlgebra 
from copy import *; 
import random; 
import timeit; 
#from ga import *;  # GAlgebra 
#from mv import *;  # static GAlgebra wrappers 
from GA_multor import *;  # ClifFred 

FGD_version = 3.1;  # update !! 

# Normalise versor (magnitude assumed positive) : 
#   version below fails to negotiate fl pt error? 
# def normalise (X) : return X/X.norm() # end def 
#   (symmetric) workaround below suffices for approximate versors! 
#def mag2 (X) : return sympify(grade((rev(GA.mv(X))*X), 0)) # end def  # coerce! 

def normalise (X) : 
  return GA.mul(X, GA.bld([1.0/sqrt(abs(GA.mag2(X)))])) # end def 


# Versor transforming orthonormal  Cl(p,q,r)  frame F to G ; 
#   L-to-R composition, frame length = n ; 
#   optional verbosity and spin continuity with previous result 
# (fails to detect mismatched isotropic  --- FIX ?? ) 
def frame_transform_versor (F, G, verb = False, Z0 = 0) : 
  # local s, t, disc, k, i, n, r, Hk, H, Z, Rk, R2, err, m1, m2, n1, n2;  
  # global GA, n, eps; 
  if verb : 
    print; print "frame_transform_versor: F, G"; print F; print; print G; print; 
    print "k, s, ||R1||, ||R2||, R, H"; print; # end if 
  
  Z = GA.bld([1.0]);  t = 1;  # main loop 
  for k in range(0, GA.dim) : 
    #Hk = GA.mul(GA.bld([t]), GA.gra(GA.mul(GA.mul(GA.rev(Z), F[k]), Z), 1)); 
    Hk = GA.mul(GA.bld([t]), GA.form(F[k], Z, 1, 1)); 
    R1 = GA.sub(G[k], Hk); m1 = GA.mag2(R1); 
    R2 = GA.add(GA.bld([GA.mag2(G[k])]), GA.mul(Hk, G[k])); m2 = GA.mag2(R2); 
    n1 = 0 if m1 == 0 else min(abs(m1), 1.0/abs(m1)); # end if 
    n2 = 0 if m2 == 0 else min(abs(m2), 1.0/abs(m2)); # end if 
    if max(n1, n2) < eps :  # copy (signatures clash or both isotropic) 
      Rk = GA.bld([1.0]); s = 1; 
    elif n1 > n2 :  # reflect 
      Rk = normalise(R1); s = -sign(m1); 
    else :  # rotate 
      Rk = normalise(R2); s = +sign(m2); # end if 
    Z = GA.mul(Z, Rk);  t = t*s;  # accumulated transform & sign 
    if verb : 
      print k+1, s, m1, m2; print; print Rk; print; print Hk; print; # end if end for 
  # finally  H = (1/Z) F Z ~ +/-G , Z = prod_k R ; 
  
  # fix frame signs if possible ( p  even,  r  zero;  p  odd,  r  nonzero);  
  #   misses mismatched isotropic generators! 			 --- FIX ?? 
  if t < 0 : Z = GA.mul(Z, GA.J);  # end if 
  
  if Z0 <> 0 :  # ensure  sign(Z)  continuous 				 --- ?? 
    Z_Z = GA.mul(Z0, GA.rev(Z));  # (old Z)/(new Z) ~ disc  (monomial) 
    disc = round(GA.lis(Z_Z)[0]); 
    if disc == 0 and n%2 == 1 : 
      disc = GA.rev(GA.dual(GA.bld([round(GA.lis(GA.dual(Z_Z))[0])]))); 
    else : disc = GA.bld([disc]); # end if 
    if GA.mag2(disc) <> 1 : 
      print "frame_transform_versor: discontinuity = ", disc; print; 
    else : Z = GA.mul(disc, Z); # end if 
  
  err = 0.0;  # roundoff error 
  H = [ GA.form(F[k], Z, 1, 1) for k in range(0, GA.dim) ]; 
  for k in range(0, GA.dim) : 
    err = err + abs(GA.mag2(GA.sub(G[k], H[k]))); # end for 
  err = sqrt(err/GA.dim); 
  if err > eps : 
    print "frame_transform_versor: error = ", err; print; # end if 
    
  if verb : 
    print "G, H, Z"; print G; print; print H; print; print Z; print; print "t, err ", t, err; print; # end if 
  return Z; # end def 

# Convert vector to versor 
def vector_to_versor (V) :  # local j, W; global n, gene; 
  W = GA.bld([0]); 
  for j in range(0, GA.dim) : 
    W = GA.add(W, GA.mul(GA.bld([V[j]]), GA.gen(j+1))); # end for 
  return W; # end def 

# Convert normalised versor to matrix: rotation angles doubled; 
#   L-to-R versors, L-to-R matrices; 
#   row  i  equals vector action of versor  X  on  e_i ; 
def versor_to_matrix (X) :  # local i, j; 
  return Matrix([GA.form(GA.gen(i+1), X, 1)[1] for i in range(0, GA.dim)]); 
  # end def 

# Convert matrix to versor: rotation angles halved; calls frame_transform_versor() ! 
#   L-to-R versors, L-to-R matrices; optional spin continuity 
def matrix_to_versor (A, Z0 = 0) :  # local i; global n, gene; 
  return frame_transform_versor([ GA.gen(i+1) for i in range(0, GA.dim) ], 
    [ vector_to_versor(A.row(i)) for i in range(0, GA.dim) ], 
    False, Z0); # end def 

# Build Givens' rotation, boost, translation in  ij-plane 
def givens(i0, j0, t, s) :  # local i, j, R, r0; global n; 
  r0 = sqrt(abs(GA.gensig[j0]*s**2 + GA.gensig[i0]*t**2));  
  if r0 < eps+eps : r0 = 1;  # both isotropic? 
  R = [ [ 1-abs(sign(i-j)) 
    for j in range(0, GA.dim) ] for i in range(0, GA.dim) ];  # unit matrix 
  R[i0][i0] = t/r0; R[j0][j0] = t/r0;  # cos(angle)  etc. 
  R[i0][j0] = -s/r0; R[j0][i0] = GA.gensig[i0]*GA.gensig[j0]*s/r0;  # sin(angle)  etc. 
  return Matrix(R); # end def 

# Factorise O(p,q,r) matrix into Givens' rotations  B = +/- R_1 ... R_m , 
#   with optional verbosity; L-to-R matrices 
def givens_factor_matrix (B, verb = False) : 
  #local A, C, R, i, j, k, r, s, t, n, m, err; global eps; 
  n = GA.dim; m = n*(n-1)/2; 
  if verb : print; print "givens_factor_matrix: B"; print B; print; # end if 
  
  C = Matrix([ [ 1-abs(sign(i-j)) 
    for j in range(0, n) ] for i in range(0, n) ]);  # unit matrix 
  A = deepcopy(B); # initially  A = B , C = 1 
  R = [ C for i in range(0, m) ];  # rotation factors 
  k = m;  # move counter 
  if verb : 
    print "k, i, j, R, A"; print; # end if 
  
  for j in range(0, n) : 
    for i in range(0, j) : 
      k = k-1;  # main loop: next move 
      R[k] = givens(i, j, A[i, i], A[i, j]); 
      A = A*R[k];  # new A_ij ~ 0 , A_jj > 0 
      R[k][i, j] = -R[k][i, j];  R[k][j, i] = -R[k][j, i]; 
      C = R[k]*C;  # update product with inverse G-rot 
      if verb : 
        print k+1, i+1, j+1; print; print R[k]; print; print A; print; 
      # end if end for end for 
  # finally  A ~ diag(1, ..., 1, +/-1) , C ~ B (except row  n ?) 
  
  err = 0.0;  # roundoff error 
  for i in range(0, n) : 
    for j in range(0, n) : 
      err = err + (B[i, j] - C[i, j])**2; # end for end for 
  err = sqrt(err/n); 
  if err > eps : print "givens_factor_matrix: error = ", err; print; # end if 
  
  if verb : print "B, R[1] ... R[m]"; print B; print; print C; print; print "err", err; print; # end if 
  return R; # end def 

# Factorise  Cl(p,q,r)  versor into Givens' rotations: L-to-R composition; 
#   calls givens_factor_matrix();  Y = R_1 ... R_m (e_1 if odd) 
def givens_factor_versor (Y, verb = False) : 
  # local k, B, Z, M, Rmul, Rmat, n, m, sig, err; global GA; 
  n = GA.dim; m = n*(n-1)/2; 
  if verb : print; print "givens_factor_versor: Y"; print Y; print; # end if 
  if GA.lis(Y)[0] <> 0 : M = GA.bld([1]);  # even grade 
  else : M = GA.gen(1);  # odd grade:  e_1  non-isotropic? 
  B = versor_to_matrix(GA.mul(Y, M)); 
  Rmat = givens_factor_matrix(B); 
  Rmul = [ matrix_to_versor(Rmat[k]) for k in range(0, m) ];  # un-pruned 
  Rmul = [ GA.add(GA.gra(Rmul[k], 0), GA.gra(Rmul[k], 2)) for k in range(0, m) ]; 

  Z = GA.rev(GA.mul(Y, M));  # 1/(Y M) R_1 ... R_m  ~  1 ,  M in {1, e_1} 
  for k in range(0, m) : 
    Z = GA.mul(Z, Rmul[k]); # end for 
  sig = GA.lis(Z)[0]; Rmul[m-1] = GA.mul(GA.bld([sig]), Rmul[m-1]);  # adjust spin 
  err = sqrt(abs(GA.mag2(GA.sub(Z, GA.bld([sig])))/2)); 
  if err > eps : print "givens_factor_versor: error = ", err; print; # end if 
  
  if verb : 
    print "R_1, ..., R_m"; print Rmul; print; 
    print "(1/Y) R_1 ... R_m, M"; print Z; print M; print; print "sig, err ", sig, err; print; # end if 
  return Rmul; # end def 

################################################################################

# Random real: component range (-1..+1) ; precision ~ 10E-16 
def rand_unit () : 
  return random.uniform(-1, +1); # end def 

# Random row vector (non-isotropic) 
def rand_vector () :  # local j, ranlis; 
  ranlis = [ rand_unit() for j in range(0, GA.dim) ]; 
  return Matrix([ [ ranlis[j] for j in range(0, GA.dim) ] ]); # end def 

# Random orthonormal  Cl(p,q,r)  product of  l  vectors, avoiding isotropic 
def rand_versor (l) :  # local i; 
  if l%2 == 0 : X = GA.bld([1]);  
  else : X = vector_to_versor(rand_vector()); # end if 
  for i in range(0, floor(l/2)) : 
    X = GA.mul(X, GA.add(GA.bld([1]), GA.mul(vector_to_versor(rand_vector()), 
      vector_to_versor(rand_vector())))); # end for 
  return normalise(X); # end def 

# Random orthonormal  O(p,q,r)  matrix product of  l  reflections 
#   Q & D : for  n  odd  &  l  odd, version below yields  det = +1 ? 
def rand_ortho (l) : 
  return versor_to_matrix(rand_versor(l)); # end def 

# Instantiate Clifford algebra, given signature list comprising  +1,-1,0's 
#   GAlgebra wrapper 						 --- UNUSED ?? 
def instantiate_GA (sigs0) :  # local j; 
  global GA, GA_J, gene, sigs, eps, n, m; 
  global e1, e2, e3, e4, e5, e6;  # export default identifiers 
  sigs = sigs0; 
  n = len(sigs0); m = n*(n-1)/2; 
  
  if n == 2 : 
    GA = Ga('e1 e2', g = sigs); 
    (e1, e2) = GA.mv(); gene = [e1, e2]; 
  elif n == 3 : 
    GA = Ga('e1 e2 e3', g = sigs); 
    (e1, e2, e3) = GA.mv(); gene = [e1, e2, e3]; 
  elif n == 4 : 
    GA = Ga('e1 e2 e3 e4', g = sigs); 
    (e1, e2, e3, e4) = GA.mv(); gene = [e1, e2, e3, e4]; 
  elif n == 5 : 
    GA = Ga('e1 e2 e3 e4 e5', g = sigs); 
    (e1, e2, e3, e4, e5) = GA.mv(); gene = [e1, e2, e3, e4, e5]; 
  elif n == 6 : 
    GA = Ga('e1 e2 e3 e4 e5 e6', g = sigs); 
    (e1, e2, e3, e4, e5, e6) = GA.mv(); gene = [e1, e2, e3, e4, e5, e6]; 
  else : 
    print; print "You're on your own, sunshine!  n = ", n; print; # end if 
  
  GA_J = 1;  # quasi-pseudar 
  for j in range(0, n) : 
    if sigs[j] <> 0 : GA_J = GA_J*gene[j]; # end if end for 
  return None; # end def 

# Verbose test suite 
def test_main (verb = True) : 
  global eps;  # local Amat, A, Bmat, B, Z, R, Y, n, m; 
  print; print "test_main(): verbosity ", verb; print; 
  n = GA.dim; m = n*(n-1)/2; 
  eps = 1.0**(n-16);  # rounding error bound (double-length floating-point) 

  # Source and target frames  A,B  via random orthonormal matrices 
  Amat = rand_ortho(floor(n/2)*2);  # right-handed to right-handed 
  A = [ vector_to_versor(Amat.row(i)) for i in range(0, GA.dim) ]; 
  Bmat = rand_ortho(floor(n/2)*2); 
  B = [ vector_to_versor(Bmat.row(i)) for i in range(0, GA.dim) ]; 
  Z = frame_transform_versor(A, B, verb); 
  Amat = rand_ortho(floor(n/2)*2);  # right-handed to left-handed 
  A = [ vector_to_versor(Amat.row(i)) for i in range(0, GA.dim) ]; 
  Bmat = rand_ortho(floor(n/2)*2+1); 
  B = [ vector_to_versor(Bmat.row(i)) for i in range(0, GA.dim) ]; 
  Z = frame_transform_versor(A, B, verb); 
  
  # Target  B  random orthonormal matrix 
  B = rand_ortho(floor(n/2)*2);  # direct 
  R = givens_factor_matrix(B, verb); 
  B = rand_ortho(floor(n/2)*2+1);  # reflective: fails! 
  R = givens_factor_matrix(B, verb); 
  
  # Target  Y  random (normalised?) versor 
  Y = rand_versor(floor(n/2)*2);  # direct 
  R = givens_factor_versor(Y, verb); 
  Y = rand_versor(floor(n/2)*2+1);  # reflective: fails for odd  n ! 
  R = givens_factor_versor(Y, verb); 
  
  return "DONE"; # end def 

# Demo with  n = 2  set globally: versor spin; 
#   versor -> matrix homomorphism; spin continuity 
def spin_disc (l = 20) : 
  # local k, pile, R, R_k, X, A, R_1, R_11, R_2, R_22; global gene; 
  print; print "spin_disc(): demos, l = ", l; print; 
  
  # ambiguous sign of versor, due to double cover of continuous rotation: 
  #   after  l  steps,  X_l ~ X_0  but  R_l ~ -R_0 ~ -1 ; 
  print "k, X_k, R_k"; print; 
  pile = float(pi/l); 
  R = GA.add(GA.bld([cos(pile)]), GA.mul(GA.bld([sin(pile)]), GA.mul(GA.gen(1), GA.gen(2))));  # rotation angle 2pi/l 
  R_k = GA.rev(R); 
  X = rand_versor(1); 
  for k in range(0, l+3) : 
    R_k = GA.mul(R_k, R); X = GA.form(X, R, 1, 1); 
    print k; print X; print R_k; print; # end for end def 
  
  # spin continuity enforced during matrix conversion to versor; 
  #   unidirectional versor product mapping to matrix product 
  A = versor_to_matrix(R); 
  R_1 = matrix_to_versor(A, GA.bld([1])); # as  k = 1 
  R_11 = matrix_to_versor(A, GA.bld([-1])); # as  k = l+1 
  R_2 = matrix_to_versor(A*A, R_1); # as  k = 2 
  R_22 = matrix_to_versor(A*A, R_11); # as  k = l+2 
  print "A, R_1, R_(l+1), R_2, R_(l+2) "; print; print A; print; print R_1; print R_11; print R_2; print R_22; print; 
  
  return "DONE"; # end def 

verbose = True; demons = True;  # print & execute switches, True / False
secs = timeit.default_timer(); 
if demons : 
  print; print "Python/ClifFred frames_givens_demo version", FGD_version; 
  GA = ClifFred([1,1,1,1]);  # 3-sphere geometry  Cl(4) 
  print test_main(verbose);  # verbose tests 
  GA = ClifFred([1,1]);  # circle geometry  Cl(2) 
  print spin_disc(20); # spin demo 
  GA = ClifFred([1,1,-1,-1,0,0]);  # mixed degenerate  Cl(2,2,2) 
  print test_main(verbose);  # verbose tests 
# end if 
secs = timeit.default_timer() - secs; 
print "Elapsed time in secs ", secs; 

################################################################################

# TODO --- 
# Wrap non-static properties: matrix.T , .row, list.append; unit matrix, etc ?? 
# spin_disc() : demo dual continuity ?? 
# test_main () : tidy identifiers ?? 
# rand_versor() : fix sign quibble? 
# frame_transform_versor() : 
#   misses mismatched isotropic generators! Maybe related to  --- INVESTIGATE ?? 
# givens_factor_versor() : 
#   Cl(2,2,2)  often reports large errors, tho matrix OK ! 	 --- BUG ?? 
# givens_factor_matrix() : 
#   Now reports large error for reflective case (formerly suppressed)! 
#   Query: test translation in 1-space CGA = GA(2, 1) ? 
#   Query: givens_factor_versor() succeeds though givens_factor_matrix() fails! 
# Wrappers & times for GA_multor; GAlgebra, Clifford ... ?? 
#   Time 18.75 sec -> 6.47 : compare GAlgebra ~141.0 / 21.8 ; printing ~0.1 sec?  
# Set 1 & unit vectors ??  Use  GA.form()  ?? 
# Document: switches, wrappers !! 
# Link switches to command line? 
# Use  NumPy  matrices, multi-prec? 
# Use mullis, addlis; from_zero() , from_prop() ?? 

