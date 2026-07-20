This is the README file for the new version of the atom simulator. The following updates are made:
1) Improved Nucleus generation using my own algorithm.
2) Improved movement physics and speed display.

ABOUT THE NUCLEUS GENERATION ALGORITHM:
Derivation: We know that the sum of the volumes of all the nucleons is going to be less than the tight-closed sphere they are in. Therefore,
we can set the following approximate equation: n * 4/3 * pi * r^3 = 4/3 * pi * R^3 where r is the radius of the nucleons, n is the number of
nucleons and R is the radius of the larger sphere. Simplifying the equation, we have: R^3 = n * r^3. Taking cube roots of both sides, we
have: R = r * n^{1/3}. Now, let us postulate that if n is odd, the centre will be occupied by a nucleon and if not, it will not be occupied.
So, R = 2xr + y where y = r if n is odd and y = 0 if n is even and x is the number of nucleons in a single line or axis. Solving for x gives
: x = (R - y) / 2r. Now, we just have to loop over the factor angles of 90 degree on the polar and azimuthal angles and use the spherical-to
-cartesian coordinate conversion formulas, we will have the unit-vector of the line on which x number of nucleons are present. To find the
coordinates of the nucleons, we just have to loop over j and evaluate 2jr + y - r = r(2j- 1) + y and multiply this by the unit vector until
we get n valid positions. Valid in the previous sentence means that none of the nucleons overlap due to similar or too close coordinates. You
might be wondering that the radius of the nucleons I draw is greater than r and so I intentionally make them overlap. This is to show respect
to the real physical laws that govern the structure of nuclei in atoms as nucleons can absolutely overlap with each other's wavefunction to
strengthen the strong nuclear attraction. I also have another rule, i.e., to draw the mirror point in the opposite quadrant just after
drawing the one in the original quadrant to preserve symmetry.

There might be flaws in my algorithm but who likes sticking to those complex, grid and other algorithms? I thought it was better writing my
own.
