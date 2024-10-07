# Two_Sex_Time_Varying_MultiState_Kin

Code to implement a combination of Hall Caswell's formal demography of kinship models: volumes II,III,IV; two-sex, multi-state kinship 
subject to time-varying demographic rates.

As inputs the model needs

1) Time-and-sex-dependent and age-specific fertility rates, categorised over distinct stage/states.

2) Time-and-sex-dependent and age-specific mortality rates, categorised over distinct stage/states.

3) Time-and-sex-dependent and age-specific probabilities of moving from one stage to another

4) Time-and-sex-dependent and age-specific redistribution of newborns (see Eq (10) in Caswell II for more details)

5) A sex-ratio (assumed static for now)

6) Specification as to which kin one wishes to analyse (if not all of Focal's network)

7) Whether the stage/demographic characteristic is the particular case of parity

8) Specification as to whether the user wishes to have accumulated numbers of kin up to age of Focal, 
   or age*stage distributions of kin given a fixed age of Focal (dist_output)
   
9) Focal's sex 

10) An integer used to specific whether the user applies abridged age-classes. I.e., nc-year age-class

11) A sequence of years which correspond to the time-varying demographic rates above (e.g., seq(1990, 2070, by = nc))


As outputs the model gives 

1) If one wishes to look at accumulated kin, and "dist_output = FALSE", then for each year, the model gives:

Kin numbers -- i.e., summed over all ages of kin -- through variable "pred_no_kin", categorised by "Stage_Kin" for each age of      Focal "Age_Focal"

2) If one wishes the distribution of kin "dist_output == TRUE", then for each year we have

kin numbers by age and stage of kin, and for each age of Focal.


# Two_Sex_Time_Varying_MultiState_Kin_

Implements the above exactly the same, however does not need to call functions from other files -- all functions inculded in one 
(very very long) script


