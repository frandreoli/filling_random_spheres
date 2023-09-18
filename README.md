# Introduction
In this repository, few algorithms are dicussed to uniformly sample hard (non-overlapping) spheres inside a cuboid. Although fully tested only for cuboid shapes, the code is designed to be potentially extendable to arbitrary shapes of the enclosing volume and and dimensions.

- The most basic algorithm (hereafter labeled as **_basic_**) is defined as follows. At each step a coordinate for  point is uniformly sampled inside the enclosing shape. Then 
