# What
Debug Data for RVH

# Install

## Compile Madness
```bash
git clone https://github.com/m-a-d-n-e-s-s/madness madness_src
mkdir madness_build
cmake -D ENABLE_MPI=OFF -D BUILD_SHARED_LIBS=ON -D ENABLE_MKL=ON -B build -S madness_src
make 
```

## Compile this program
```bash
cmake .
make
```

## Make Gridpoints
```bash
python make_gridpoints.py
```

# Run
```
h2j input=input
```
will use [input](input) file provided here  
will use [gridpoints.txt](gridpoints.txt) file provided here  
will print out "[x,y,z] : J(x,y,z)  
where J is \int 1/(r-r') O(r')^2 dr' and O is the left STO-3G function
