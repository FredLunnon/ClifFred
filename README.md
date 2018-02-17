# ClifFred

  Module  <GA_multor.py>  comprises Python source for a lightweight, portable, 
reasonably efficient geometric algebra engine, implementing basic operations 
on multors (multivectors) in a degenerate Clifford algebra  Cl(p,q,r) , 
aka  |R^p,q,r ; some documentation is provided in  <GA_multor.txt> .  
Components are notionally approximate numerical (double-length floating-point), 
though with appropriate restrictions they may belong to any supported numerical 
type, such as integer, complex, or rational, symbolic (using <SymPy>), etc.  
This release should be regarded as beta-stable. 

  Module  <GA_scripts.py>  contains several examples of complete though 
sparsely documented applications using the engine.  The first section includes 
a short toolkit for generating random test data, followed by increasingly 
elaborate sections demonstrating various GA identities.  
The final two sections comprise nontrivial geometric applications: measurement 
of a tetrahedron in Euclidean 3-space, and the Apollonian problem in contact 
(Poincaré) 2-space. 
 
A more elaborate and developed application  <GA_givens.py>  documented in  
<GA_givens.txt>  implements frame transformation and Givens' decomposition. 

  Version 1.2 introduces some extra functions, tinkers with specifications 
which may necessitate occasional changes to existing programs, and attempts 
to facilitate comparisons between competing packages via wrappers. 
It is intended to provide wrappers for further packages in future; 
although meaningful comparisons are hindered by their divergent design 
targets and specification of subsidiary operators. 

