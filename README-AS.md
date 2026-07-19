Hi whoever is reading this. This is an atom simulator which can simulate atoms. The rough structure of my code is mentioned below:
1) The user enters the atomic number (1-30) which is then used to find the electronic configuration in the form of numbers.
2) Each orbital is assigned a numerical value (s -> 1, p -> 2, d -> 3)
3) Until pointcnt number of electron superpositions are received, the loop is not ended.
4) Then it calculates the probabilities of the electrons for different orbitals using predefined formulas and hardcoded
normalization values.
5) A random rejection-sampling technique similar to Monte Carlo sampling is used to evaluate whether to accept the electron superposition or
not.
6) Using raylib's mesh and model things, the millions of positions and colors are loaded on the VRAM.
7) A blue sphere is drawn at the origin to indicate the nucleus and the electron superpositions are also drawn.
8) Other side-features like FPS, crosshair, WASD navigation are also included.

NOTE: Anyone who is comprehending my code must known basic vector things and some quantum mechanics.
