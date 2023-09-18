# Introduction
In this repository, few algorithms are dicussed to uniformly sample hard (non-overlapping) spheres inside a cuboid. Although fully tested only for cuboid shapes, the code is designed to be potentially extendable to arbitrary shapes of the enclosing volume and and dimensions. Let us suppose that the hard spheres are $N$ have radius $R$. Let us suppose as well to be working in $D$ dimensions (the code was fully tested for $D\leq 3$ dimensions). The metric tested is euclidian, i.e. $|\mathbf{r}\_i-\mathbf{r}\_j|=\sqrt{ \sum\_{x=1,\dots, D} (r\_j^x-r\_j^x)^2}$. 

- The most basic algorithm (hereafter labeled as **_basic_**) is defined as follows. At each step $1\leq i\leq N$ a set of coordinates $\mathbf{r}\_i$ for the _i_-th point is uniformly sampled inside the chosen enclosing shape (the standard shape tested here is a cuboid). If the distance condition $|\mathbf{r}\_j-\mathbf{r}\_i|\leq 2R$ is satisfied for all the previous $j=1, 2, \dots, i-1$ points, then the new point is accepted, otherwise rejected.

<!--- - Second, we tested another possible approach (hereafter labeled as **_joint_**), where a set of $N$ coordinates is directly sampled from the beginning. Then,--->

- A more efficient algorithm (hereafter labeled as **_grid_**) is developed as follows. It is based on the idea of dividing the enclosing volume into small cubes of size $L$, to form a $D$-dimensional grid. For a generic case (such as a generic cuboid), we encapsulate the enclosing shape into a larger grid, and then at each step we reject the new point if this ends up laying outside the target shape. The grid constant $L=2R\sqrt{D}$ is chosen in such a way that two distinct points (i.e. the centers of two hard spheres) cannot be inside the same cube, as they would necessarily overlap. At each step _i_ of the code, the algorithm proposes a new point by first randomly choosing one among the (empty) small cubes, and then uniformly sampling a countinuous coordinate $\mathbf{r}_i$ inside. To accept the point, the algorithm checks its distance from only the points inside the neighbouring cubes, rather than from all the other points. To this aim, the cubes are considered neighbouring up to the $m$-th neighbouring order, where `m=ceil(Int, sqrt(D))`.

- 


# Scaling examples

# Kolmorogov-Smirnov tests
