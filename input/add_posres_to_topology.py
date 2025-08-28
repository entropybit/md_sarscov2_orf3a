import argparse
import numpy as np


if __name__ == "__main__":


    parser = argparse.ArgumentParser(description='Rename HOH to SOL')
    #parser.add_argument('-f','--foo', help='Description for foo argument', required=True)
    #parser.add_argument('-b','--bar', help='Description for bar argument', required=True)
    parser.add_argument("INTOP", help=" input topol.top ")
    parser.add_argument("MOLNAME", help=" molecule name for which restraint need to be added")
    parser.add_argument("DPOSRES", help=" name of -DPOSRES define used for restraint section")
    parser.add_argument("RESITP", help=" path to posre_X.itp containing restraints")
    parser.add_argument("OUTTOP", help=" output topol.top ")
    args = parser.parse_args()


    top_lines = None
    with open(args.INTOP) as f:
        top_lines = f.readlines()    

    top_lines = np.array(top_lines)

    ind_moltype=np.where(top_lines == "[ moleculetype ]\n")[0]
    ind_molnames = ind_moltype+2
    print(ind_moltype)
    print(top_lines[ind_molnames])


    check_match = [args.MOLNAME in x for x in top_lines[ind_molnames]]
    ind_idxEntry_to_add = np.where(check_match)[0][0] + 1 % len(check_match)
    
    print(check_match)
    print(ind_idxEntry_to_add)

    ind_add_here = ind_moltype[ind_idxEntry_to_add]-1


    posres_definition = np.array(
        [
            f"\n",
            f"; Include Position restraint file\n",
            f"#ifdef {args.DPOSRES}\n",
            f"#include \"{args.RESITP}\"\n",
            f"#endif\n",
            f"\n"
        ]
    )

    print(" going to add the following ::")
    print(posres_definition)
    print("")

    top_lines = np.concatenate(
        [
            top_lines[:ind_add_here],
            posres_definition,
            top_lines[ind_add_here:]

        ]
    )
    
    with open(args.OUTTOP, "w") as f_out_top:
        for l in top_lines:
            f_out_top.write(l)


