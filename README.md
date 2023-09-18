# Introduction
In this repository, few algorithms are dicussed to uniformly sample hard (non-overlapping) spheres inside a cuboid. Although fully tested only for cuboid shapes, the code is designed to be potentially extendable to arbitrary shapes of the enclosing volume and and dimensions. Let us suppose that the hard spheres are $N$ have radius $R$. Let us suppose as well to be working in $D$ dimensions (the code was fully tested for $D\leq 3$ dimensions). The metrix is euclidian, i.e. $|\mathbf{r}\_i-\mathbf{r}\_j|=\sqrt{ \sum\_{x=1,\dots, D} (r\_j^x-r\_j^x)^2}$. 

- The most basic algorithm (hereafter labeled as **_basic_**) is defined as follows. At each step $1\leq i\leq N$ a set of coordinates $\mathbf{r}\_i$ for the _i_-th point is uniformly sampled inside the chosen enclosing shape (the standard shape tested here is a cuboid). If the distance condition $|\mathbf{r}\_j-\mathbf{r}\_i|\leq 2R$ is satisfied for all the previous $j=1, 2, \dots, i-1$ points, then the new point is accepted, otherwise rejected.

<!--- - Second, we tested another possible approach (hereafter labeled as **_joint_**), where a set of $N$ coordinates is directly sampled from the beginning. Then,--->

- A more efficient algorithm (hereafter labeled as **_grid_**) is developed as follows. It is beased on the idea of dividing the enclosing volume into small cubes of size $L$, to form a $D$-dimensional grid. For a generic enclosing shape (such as a generic cuboid), encapsulate this shape into a larger grid, and then at each step reject the points if they end up laying outside this target shape. The grid constant $L=2R\sqrt{D}$ is chosen in such a way that two distinct points (i.e. the centers of two hard spheres) are forbidden to be inside the same cube, as they would necessarily overlap. 


# Scaling examples

# Kolmorogov-Smirnov tests
