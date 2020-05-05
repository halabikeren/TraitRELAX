############################################ How to run ############################################

# The runs need to be carried out from the 'Examples' directory. Define:

EXAMPLES_DIR="/your/path/to/Examples"
TRAIT_RELAX_EXEC="/your/path/to/TraitRELAX/traitrelax"

# And run (change '1' to '2' for the second example):

cd "${EXAMPLES_DIR}"
"${TRAIT_RELAX_EXEC}" param="${EXAMPLES_DIR}"/example1.bpp > "${EXAMPLES_DIR}"/res1_full_log.txt

###################################### Info about the examples ######################################

### Example 1

The input in example1.bpp is a replicate from the TraitRELAX simulation study.
The parameter values were:

mu=8
pi0=0.5
kappa=2
omega0=0.1
omega1=0.8
omega2=2
p0=0.5
p1=0.4
k=0.2

TraitRELAX detects a significant (p-value=0.0336) trait-related relaxation with parameter estimates:

mu=3.25541190565717
pi0=0.05
kappa=1.87347952238339
omega0=0.00490177859
omega1=0.490177859225133
omega2=1.20269784121443
p0=0.0380727606300336
p1=0.961927239
k=0.258846591077448

When using 16 cores, the estimated runtime is 2.5 hours.
You can compare your output files to those provided in "${EXAMPLES_DIR}"/data/example1/

### Example 2 (will be added soon)

This is a toy example of a super small data set (4 species and 1 codon position), designed to ensure
a valid execution of the program. It should run very quickly.
You can compare your output files to those provided in "${EXAMPLES_DIR}"/data/example2/
