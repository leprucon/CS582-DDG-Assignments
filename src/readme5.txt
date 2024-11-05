1.1: faced a bunch of bugs with compilation, but otherwise it was straightforward. I just followed documentation, and commented out lines when they weren't needed to generate each data series. I also made sure that the size of each successive matrix was evenly increased
1.2.2: pretty much the same as 1.1, but I used a random number generator to randomly pick positions. there are two pairs of elements for each random number, so I generated a new number for 5% of the size of the matrix for 10% total.
2.1: I implemented deflate by directly following the equation, and did the same for power iteration
2.2: At first I generated nice straight logarithmic graphs, but then I discovered this was due to a bug. Now my result simply does not converge, despite my best efforts to proof read my logic.
3.1: I derived the functions analytically before implementing them in code. 
3.2: I decided to apply the same trick as before, and use decomposition to more quickly calculate the inverse hessian term.
