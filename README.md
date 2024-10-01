# Two_Sex_Time_Varying_MultiState_Kinship

Code to implement a combination of Hall Caswell's formal demography of kinship models: volumes II,III,IV; two-sex, multi-state kinship 
subject to time-varying demographic rates.

As data inputs the model needs, separate for males and females

1) Time-dependent and age-specific fertility rates, categorised over distinct stage/states.

2) Time-dependent and age-specific mortality rates, categorised over distinct stage/states.

3) Time-dependent and age-specific probabilities of moving from one stage to another

4) A sex-ratio of male:female newborns (alpha)

For now, male and female treated identical for

5) Redistributing newborns from each age*stage combination to a new age*stage

6) A sequence of years corresponding to the timescale of the demographic rates 

As theoretical choices, the model needs

1) Choice as to whether Focal is male or female (default is female)

2) Choice as to whether the specific state/stage is parity

3) Choice (if not parity) of the stage in which Focal is born


An example of the model is given using UK parity as a state
