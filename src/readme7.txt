Question 1 (For README): Iterate your method on the three provided cylinders. One
of them shows a behavior that is different from the other two. Can you explain what is
happening? Is the result consistent with the goal of the minimal surface optimization?
 - The reason why this cylinder 'implodes' is because the cylinder is too tall- conceptualizing it as two rings with a soap bubble between them, drawing the rings farther and farther apart will eventually cause the walls of the bubble to collapse in on itself. 
 - The result isn't consistent with minimal surface optimization, because it begins taking erratic values after a couple iterations, causing the mesh to not converge, and eventually causing a segfault.

Question 2 (For README): Does the same effect happen if you replace the cotan Laplacian with the uniform Laplacian? Elaborate your answer
 - The same effect no longer occurs when you use the uniform laplacian.
 - The uniform laplacian also no longer adheres to the minimal surface calculation