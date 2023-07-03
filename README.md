# uGBE_BCC 

Universal grain boundary function for BCC metal was estimated via based on RSW function.
The function is computed though MATLAB code version R2022a.

To estimate grain boundary energy, the inputs consist of 
1. input data file (fnData) that contain 18 columns of the orientations in the two grains (P and Q) 
2. coherent twin grain boundary energy (GBE_coh_twin) in unit of J/m<sup> 2 </sup>
3. type of 'alkali' or 'transition' BCC metals (type).

Command 
```
uGBE(fnData, GBE_coh_twin, type)
```

# Example
We have example file in this repository ("Example_GB.xlsx").
This file is information of the sigma 23 (310) boundary, it has 18 columns of the orientation in the two grains (P and Q).

     |                  Matrix P                  |                  Matrix Q                  |
     | x1 |    |    | y1 |    |    | z1 |    |    | x2 |    |    | y2 |    |    | z2 |    |    |
     |  6 |  2 |  0 |  3 | -9 |  5 |  2 | -6 |-12 |  6 |  2 |  0 |  3 | -9 | -5 | -2 |  6 |-12 |

Our previous work has reported grain boundary energy in Fe which has coherent twin grain boundary energy of 0.262260282 J/m<sup> 2 </sup>, and type of 'transition'.

```
uGBE('Example_GB.xlsx', 0.262260282,'transition')
```

Function returns 1.2920 as an energy value of sigma 23 (310) boundary in unit of J/m<sup> 2 </sup>.
