Coding:
  Computing Vertex Normals:
    - Constant weights was straightforward. I used face normals given by geometry central
    - Weighted by area was also straightforward. I used face areas given by geometry central
    - For weighting by incident angles, I directly used the vertex normals computed by geometry central, because I read in the documentation that they computed vertex normals with incident angle weights.
  2.1:
    - I essentially followed the discretized operator equation to a T, except I decided to set boundary vertex weights to zero arbitrarily. What I was supposed to do with boundary vertices was unclear.
  2.2:
    - I ran into a small hiccup with figuring out how to associate adjacent edges with adjacent vertices. While my solution works, the logic feels clunky.
  2.3:
    - Was straightforward, just followed the formula.
  

Small vs. Larger sphere:
  Visually, for both the Uniform Laplacian and Laplace-Beltrami operators, the small and larger spheres look the same- the resulting texture is just scaled up on the larger sphere proportionally (which is expected, since the color range is scaled by max/min values).
  When we observe the min/max curvature values however, we see that the Uniform Laplacian has all of its curvatures scaled down by a factor of 10 from the large to small sphere. For the Laplace-Beltrami operator, this is reversed- the curvatures are scaled up by a factor of 10 from the large sphere to the small sphere.
