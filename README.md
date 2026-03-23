Repository for our solutions of CGSHOP2026: Central Triangulation under Parallel Flip Operations.
 

For further information, see the paper "ETH Flippers Approach to Parallel Reconfiguration of Triangulations: SAT formulation and Heuristics" (L. Battini, M. Milenković, 2026) on arxiv (to be released).

See the XOR-SAT folder for the XOR-SAT formulation, and how to use it. 
The two\_triangulations folder contains an implementation of the paper "Simultaneous edge flipping in triangulations" (J Galtier, F Hurtado, M Noy, S Pérennes, J Urrutia, 2003), which was not used in the competition because the resulting solutions are not good.

Dependencies
- Qt 6 for visualization (it should compile without it and just visualize nothing)
- CGAL 6.0.1

<h2>Build</h2>

```
cmake .
```

```
make
```

<h2>Run</h2>

```
./centralTriangulation pathToInputFile pathToOutputFile methodName verbose/nonVerbose
```
where method refers to how the solution is computed:
- generateCNFNU requires additional parameter pathToidxFile listOfDistancesToTheCenter (after all the other parameters). idxFile stores information which is used to convert the assignment file to solution file (see next command), listOfDistancesToTheCenter is a series of integer, on per input triangulation, in the same order as in the input file. This method generates a SAT formulation in the DIMACS format, which is satisfiable if and only if it there is a flip sequence transforming each input triangulation into a common center with the number of flips given by the corresponing distance to the center parameter.
- assignmentToFlipListNU, converts an assignment file (produced by a SAT solver) to solution in the competition .solution.json format. It requires additional parameter idxFile (at the end).
- flipDiff, which is the simple heuristic of only flipping 
diagonals that are not shared between two triangulations (this only works for very simple inputs);
- reduceIts, a heuristic which flips diagonals when that reduces the number 
of intersections
- reduceItsCenter, same as reduceIts, but it tries more possible centers
- by default, if no method is specified, it uses Delaunay as a central triangulation. This method is not kept 
up to date (it doesn't use the abstract class) because initial tests showed some terrible results
(in particular, it is not possible to "just add" parallel flips to the Lawson flip algorithm)


<h2>Scripts</h2>
multipleInputs.py computes an heuristic solution for multiple files.
searchSAT\*.py are various scripts that compute exact solution, by searching over the solution space, and verifying each instance with generateCNFNU, followed by a call to a SAT solver, and finally a call to assignmentToFlipListNU.
verifyAndVisualize.py uses the verifier provided by the competition organizers, see [github.com/CG-SHOP/pyutils26](https://github.com/CG-SHOP/pyutils26), to verify if a solution is valid, print the objective and optionally visualize the solution.
Various parameters for those scripts (such as the path to the SAT solver) need to be specified in the source code.

<h2>Adding another solution method</h2>
flippingDS is an abstract class that should already provide all the methods needed to
implement various solution methods. 

geometryUtils provides some additional helper methods, eg the one that returns whether
a diagonal can be flipped.

inOut is used to parse inputs, and write output files as json. 
Currently, it returns a CGAL constrained triangulation as 
well as other useful data structures, and these can then be used to initialize flippingDS.

See differentDiagonals.cpp and main.cpp for an example usage.
Remember to update CMakeLists.txt as well when adding a new file.
