1.2:
* Do you experience an analogous behavior for surfaces? Experiment with various meshes,
time steps, and number of iterations. Briefly comment on your observations (no formal
proof expected).
 - (using 0.5 timestep) I found that for surfaces, applying smoothing operations for even 5000 iterations would not bring the surface to a convex shape. Extremities of surfaces (such as the ears of the bunny model) would converge to a sharp point with concave curves. Though eventually the surfaces grow smaller and smaller, converging to a point.

1.3:
* Why are the methods unstable when applied on a bad quality mesh?
 - As mentioned in lecture, the uniform laplace-beltrami discretization is a "bad approximation for irregular triangulations".
 - We can see that the Bad Max has many irregular triangulations, while the Nice Max has regular triangulations (closer to equilateral)

* In your submission, please also provide the following:
• At least three images of interesting feature enhancement results.
• In the readme.txt file, explain how you achieve these results (including parameter values, the order of applying enhancement/smoothing, the iterations of enhancement/smoothing, etc.).
 - For all 3 feature enhancement images, I used an F.E. coefficient of 10. I first apply one or two iterations of uniform laplacian smoothing to resolve some mesh triangulation irregularities. 
 - Then I apply feature enhancement with cotangent weights until I visually deem the mesh to have become too 'spiky' or with overlapping triangles, which I ameliorate by running 1-5 iterations of uniform smoothing until all 'spikes' are removed. 
 - I repeatedly alternate these feature enhancement and smoothing steps until I get an interesting result.

* From a signal processing point of view, what are the effects of the operations you
apply, and why do they produce the results you obtain? Please provide your answer in the readme.txt file.
 - From a signal processing PoV, we are effectively increasing the magnitude of high-frequency details to make them more pronounced, while keeping the lower-frequency components unaffected.

2.0:
* Replace the squared Discrete Laplace-Beltrami matrix L^2 by the Discrete Laplace-Beltrami
matrix L (Figure 8). Comment on the difference between the results.
 - I tried both deformation equations on the sphere model, and found that L^2 keeps the sphere 'elliptical' or convex when stretched, while L causes the fixed/displaced points to by pulled out in a spike, leading to concave surfaces.

* Replace the cotangent weights used in Discrete Laplace-Beltrami matrix L with uniform
weights. Comment on the difference between the results.
 - The only difference I noticed (when testing on the plane mesh) was that the deformation with uniform weights causes the deformation to become a lot less smooth.

* Comment on how the computed deformations differ from a real physical material, e.g.an
elastic rubber membrane or a thin metal sheet. Create your own example shapes and deformations where those differences would be easily explainable. Please upload your(.obj/.stl)
model as well as images of its deformations in the submission.
 - Unlike real physical materials, the computed deformations don't preserve mass, can cause the surface to clip through itself, and can stretch/bend indefinitely without compromising the surface 'material'. 
 - There are also certain edge cases where arrangements of fixed points can cause the surface to behave strangely after deformation. For example, on a sphere, having a fixed point in the middle of a ring of fixed points will cause the surface to form a dimple at the middle fixed point when displaced outwards on the other side. When a single point on the sphere is displaced outwards, the opposite side of the sphere is displaced in the opposite direction.
