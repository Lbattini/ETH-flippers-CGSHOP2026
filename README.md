Repository for our solutions of CGSHOP2026: Central Triangulation under Parallel Flip Operations 

For further information, see notes/notes.pdf

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
./centralTriangulationClion pathToInput pathToOutput methodName
```
where method refers to how the solution is computed. Currently, there are two methods
- flipDiff, which is the simple heuristic of only flipping 
diagonals that are not shared between two triangulations (only works for very simple inputs);
- reduceIts, a heuristic which flips diagonals when that reduces the number 
of intersections
- reduceItsCenter, same as reduceIts, but try more possible centers
- by default, if no method is specified, it uses Delaunay as a central triangulation. This method is not kept 
up to date (it still doesn't use the abstract class) because initial tests showed some terrible results
(in particular, it is not possible to "just add" parallel flips to Lawson flip algorithm)


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
