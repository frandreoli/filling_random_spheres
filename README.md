j# Introduction
In this repository, few algorithms are dicussed to uniformly sample hard (non-overlapping) spheres inside a cuboid. Although fully tested only for cuboid shapes, the code is designed to be potentially extendable to arbitrary shapes of the enclosing volume and and dimensions. Let us suppose that the hard spheres have radius R. 

- The most basic algorithm (hereafter labeled as **_basic_**) is defined as follows. At each step _i_ a set of coordinates $\mathbf{r}_i$ for the _i_-th point is uniformly sampled inside the chosen enclosing shape (the standard shape tested here is a cuboid). If the distance condition $|\mathbf{r}_j-\mathbf{r}_i|\leq 2R$ is satisfied for all the previous $j=1...i-1$ points, then the new point is accepted, otherwise rejected.

- Ano
