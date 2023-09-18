# Introduction
In this repository, a grid-based algorithm is developed to uniformly sample hard (non-overlapping) spheres inside a cuboid. Although fully tested only for cuboid shapes, the code is designed to be potentially extendable to arbitrary shapes of the enclosing volume and dimensions. 

Hereafter, we suppose that the hard spheres are $N\_{\text{spheres}}$ have radius $R\_{\text{spheres}}$, and that we are working in $D$ dimensions (the code was fully tested for $D\leq 3$ dimensions). The chosen metric is euclidian, i.e. $|\mathbf{r}\_i-\mathbf{r}\_j|=\sqrt{ \sum\_{x=1,\dots, D} (r\_j^x-r\_j^x)^2}$. 

- **Basic algorithm.** The most basic algorithm (hereafter labeled as **_basic_**) is defined as follows. At each step $1\leq i\leq N$ a set of coordinates $\mathbf{r}\_i$ for the _i_-th point is uniformly sampled inside the chosen enclosing shape (the standard shape tested here is a cuboid). If the distance condition $|\mathbf{r}\_j-\mathbf{r}\_i|\geq 2R\_{\text{spheres}}$ is satisfied for all the previous $j=1, 2, \dots, i-1$ points, then the new point is accepted, otherwise rejected.

<!--- - Second, we tested another possible approach (hereafter labeled as **_joint_**), where a set of $N\_{\text{spheres}}$ coordinates is directly sampled from the beginning. Then,--->

- **Grid-based (exact) algorithm.** A more efficient algorithm (hereafter labeled as **_grid_**) is developed as follows. It is based on the idea of dividing the enclosing volume into small cubes of size $L\_{\text{grid}}$, to form a $D$-dimensional grid. For a generic case (such as a generic cuboid), we encapsulate the enclosing shape into a larger grid, and then at each step we reject the new point if this lays outside the target shape. The grid constant $L\_{\text{grid}}=2R\_{\text{spheres}}\sqrt{D}$ is chosen in such a way that two distinct points (i.e. the centers of two hard spheres) cannot be inside the same cube, as they would necessarily overlap. At each step _i_ of the code, the algorithm proposes a new point by first randomly choosing one among the (empty) cells of the grid, and then uniformly sampling a countinuous coordinate $\mathbf{r}_i$ inside this latter. To accept the point, the algorithm checks its distance $|\mathbf{r}\_j-\mathbf{r}\_i|\geq 2R\_{\text{spheres}}$ from the points $j$ that occupy the neighbouring cubes, rather than from all the accepted points. To this aim, the cubes are considered neighbouring up to the $m$-th neighbouring order, where `m=ceil(Int, sqrt(D))`.

- **Grid-based (approximated) algorithm.** The previous algorithm can be made more efficient by paying the price of introducing approximations. Specifically, the grid size can be fixed to $L\_{\text{grid}}=2R\_{\text{spheres}}$, while still assuming that each small cube can only contain one point. The error probability can be estimated as the probability that two points _i$ and _j_ randomly sampled inside such cube of size $L\_{\text{grid}}=2R\_{\text{spheres}}$ have an acceptable distance $|\mathbf{r}\_j-\mathbf{r}\_i|\geq 2R$. Numerically, we found this error probability to be $\sim 3$\% in 2D and $\lesssim 10$\% in 3D, while being zero in 1D. This choice of $L\_{\text{grid}}=2R\_{\text{spheres}}$ allows to speed up the code by constructing less grid cubes and, at the same time, fixing the neighbouring order to $m=1$.


# Scaling example
Here, we provide an example of the performance of the three algorithms in 3D. To this aim, we fix the enclosing shape to be a cubic box of sizes $L_{\text{tot}}=1$ and we define the filling fraction $f$ (or packing density) as 

$$f= \dfrac{ V\_{\text{spheres}}N}{ V\_{\text{tot}}} =  \dfrac{ 4\pi R\_{\text{spheres}}^3N}{ 3 L_{\text{tot}}^3}. $$

According to the Kepler conjecture (see [Sphere Packing Problem](https://mathworld.wolfram.com/SpherePacking.html)), the maximum packing density in 3D is $f\_{\text{max}} = \pi/(3\sqrt{2})\simeq 0.74$. 

Here below we show the scaling results when fixing $f=0.4f\_{\text{max}}$. We tested the CPU time spent on the code by averaging over $\sim 20000$ repetitions, and spanning values of $10\lesssim N\_{\text{spheres}}\lesssim 1000$. We found the empirical scalings:

- _Basic_ algorithm: $\langle T \rangle \sim N\_{\text{spheres}}^{2.3}$.
- _Grid_ algorithm: $\langle T \rangle \sim N\_{\text{spheres}}^{1.8}$.
- _Grid, approximated_ algorithm: $\langle T \rangle \sim N\_{\text{spheres}}^{1.5}$.


<p align="center">
 <img width="500" height="333" src="https://github.com/frandreoli/filling_random_spheres/assets/37184096/5333ca50-a968-4be8-9c0b-171d67a3cfe9">
</p>

  


<!--- 
---
![git_0 4_test](https://github.com/frandreoli/filling_random_spheres/assets/37184096/5333ca50-a968-4be8-9c0b-171d67a3cfe9) 
<p align="center">
  <em>Fig. 1: Scaling of the three algorithms, given f=0.4.</em>
 <span class="math display">\[y = \frac{a}{b} + c^2 + d\]</span>
</p>
--->

 

# Kolmorogov-Smirnov tests
