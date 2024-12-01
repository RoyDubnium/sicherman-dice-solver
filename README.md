# sicherman-dice-solver
A program that generates all valid combinations of dice which produce the same probability distribution as a pair of normal dice. Also included are my precomputed solutions for up to 360 sides.
The current rust program contains 3 versions, the fastest version, a low memory version, and an ultra low memory version. This is to prevent out-of-memory errors when dealing with larger dice, though keep in mind the low memory versions are slower.
