#ifndef ENSEMBLE_RES_COMP_H
#define ENSEMBLE_RES_COMP_H

#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "macros.h"
#include "smalloc.h"
#include "svm.h"
#include "tpxio.h"
#ifdef GRO_V5
#include "atoms.h"
#include "fatalerror.h"
#include "pargs.h"
#include "topology.h"
#include "trxio.h"
#else
#include "gmx_fatal.h"
#include "statutil.h"
#endif

#define LABEL1 -1 // classification label for trajectory 1
#define LABEL2 1 // classification label for trajectory 2
#define GAMMA 0.4 // default gamma parameter for svm_train
#define COST 100.0 // default C parameter for svm_train

/* Indices of filenames */
enum {eTRAJ1, eRES1, eNUMFILES};

/** Struct for holding eta data */
typedef struct {
    // Input parameters
    // See ensemble_comp function for how they are used.
    const char *fnames[eNUMFILES];
    real gamma;
    real c;
    int nthreads;
    output_env_t oenv;

    // eta output for atoms
    int natoms; // number of atoms
    atom_id *atom_IDs; // array of atom IDs, taken from fnames[eNDX1]. size = natoms

    // eta output for residues
    int nres; // number of residues
    int *res_IDs; // array of residue IDs. size = nres
    const char **res_names; // names of the residues. array size = nres
    int *res_natoms; // number of atoms per residue. array size = nres
    real *eta; // eta value of each residue. array size = nres

    // the following values may or may not be set, and are used
    // internally by ensemble_comp.

    // number of total atoms in each system
    // including the ones excluded by index files.
    int natoms_all;

    int DOT;
    int STD;
    
} eta_res_dat_t;


void init_eta_dat(eta_res_dat_t *eta_dat);
/* Initializes an eta_res_dat_t struct, such as setting pointers to NULL and setting default parameters.
 */

void free_eta_dat(eta_res_dat_t *eta_dat);
/* Frees the dynaimc memory in an eta_res_dat_t struct.
 */
 
void calc_correlations(eta_res_dat_t *eta_dat);
/* Calculates the centers of mass
 */

void atom_printer(eta_res_dat_t *eta_dat,
    rvec **x1,
    int nframes,
    t_atoms *atoms);

void calc_corr_matrix(eta_res_dat_t *eta_dat,
    rvec **x1,
    int nframes,
    t_atoms *atoms);

#endif // ENSEMBLE_RES_COMP_H