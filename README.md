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
We have example file in this repository ("Example_coherent_twin_GB.xlsx").
In this file that is information of coherent twin boundary, it has 18 columns of the orientation in the two grains (P and Q).

     |                  Matrix P                  |                  Matrix Q                  |
     | x1 |    |    | y1 |    |    | z1 |    |    | x2 |    |    | y2 |    |    | z2 |    |    |
     |  4 |  2 |  2 |  1 | -1 | -1 |  0 |  2 | -2 |  4 |  2 |  2 | -1 |  1 |  1 |  0 | -2 |  2 |

Our previous work has reported grain boundary energy in Fe which has coherent twin grain boundary energy of 0.262260282 J/m<sup> 2 </sup>, and type of 'transition'.

```
uGBE('Example_coherent_twin_GB.xlsx', 0.262260282,'transition')
```

Function returns 0.2204 as a energy value of this coherent twin boundary in unit of J/m<sup> 2 </sup>.
