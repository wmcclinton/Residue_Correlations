#include "gkut_log.h"
#include "ensemble_res_comp.h"

int main(int argc, char *argv[]) {
    const char *desc[] = {
        "res_correlation calualtes the correlation between the residues."
    };

    eta_res_dat_t eta_res_dat;
    init_eta_dat(&eta_res_dat);

    gk_init_log("eta.log", argc, argv);

    t_filenm fnm[] = {
        {efTRX, "-f1", "traj1.xtc", ffREAD},
        {efSTX, "-res", "res.pdb", ffREAD}, // provides residue information
    };

    t_pargs pa[] = {
        {"-d", FALSE, etINT, {&eta_res_dat.DOT}, "set correlation to use Dot Product"},
        {"-s", FALSE, etINT, {&eta_res_dat.STD}, "set correlation to be normalized"}
    };

    parse_common_args(&argc, argv, 0, eNUMFILES, fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &eta_res_dat.oenv);

    eta_res_dat.fnames[eTRAJ1] = opt2fn("-f1", eNUMFILES, fnm);
    eta_res_dat.fnames[eRES1] = opt2fn_null("-res", eNUMFILES, fnm);

    calc_correlations(&eta_res_dat);
    free_eta_dat(&eta_res_dat);

    gk_print_log("%s completed successfully.\n", argv[0]);
    gk_close_log();

    return 0;
}
