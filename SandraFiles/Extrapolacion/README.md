The following codes are updated after finding the errors of missing 2 in Vsb and also missing -2 in the Vhf i Vhf2 matrices. Also updated after having done the crosscheck of the way of finding the A and B errors.

## Code for the hyperfine spectrum
### Explanation of the files
- `CrosscheckIntervals.m`: This code uses the theorem of Hellmann-Feynmann to check if the errors of A and B are well computed. For each state we should be able to obtain $a,b$ and they should be the same as the ones in `IntervalsConfian.m` however they are not. We have found that linearizating was wrong!
 
$$\hat{H}|\Psi\rangle = E |\Psi\rangle \hspace{0.5cm} \rightarrow \hspace{0.5cm} E=\frac{\langle\Psi|\hat{H}|\Psi\rangle}{\langle\Psi|\Psi\rangle} = c + aA + bB \hspace{0.5cm} \rightarrow \hspace{0.5cm} \text{Theorem:} \hspace{0.5cm} \frac{\partial E}{\partial \lambda}=\frac{\langle\Psi|\frac{\partial\hat{H}}{\partial \lambda}|\Psi\rangle}{\langle\Psi|\Psi\rangle}=a,b$$

- `IntervalsConfian.m`: **OUTDATED**. This code used linealization to obtain the errors of A and B and also of the spectrum. However, after crosscheking them using `CorsscheckIntervals.m`
- `IntervalsConfian_extens.m`: **OUTDATED**
- `SaveWaveFunc.m`: Simply writes the wave functions in .txt files
- `Vhf_Ct_1.m`, `Vhf_Ct_blau.m`, `Vhf_Ct_roig.m`: **OUTDATED**. They were useful for just computing fewer states and using constants instead of $V_{hf}$ and $V_{hf2}$.
- `Vhf_Vhf2_TrobarK_blaus.m`: This is the code with which we obtain the optimal A and B for the lattice data. This is also now used to obtain their errors. We just substitute the lattice data with their value $\pm$ lattice errors to obtain $A_\pm$ and $B_\pm$. Hence, A and B errors will be `max(Aplus - A, A - Aminus)`.
- `Vhf_Vhf2_TrobarK_rojos.m`: **OUTDATED**
- `wf_plot.m`: Plots all the components of the wave functions for a determinated state. All the components mean all the $P^{LJ}_{1\mathcal{J}\mathcal{M}}$ rellevant for each matrix box in equations (8),(9),(10),(11).
