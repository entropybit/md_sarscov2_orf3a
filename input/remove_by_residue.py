import argparse


if __name__ == "__main__":


    parser = argparse.ArgumentParser(description='Removes all lines for a give residue RES from gro file')
    #parser.add_argument('-f','--foo', help='Description for foo argument', required=True)
    #parser.add_argument('-b','--bar', help='Description for bar argument', required=True)
    parser.add_argument("INGRO", help=" input GROFILE.gro")
    parser.add_argument("INTOP", help=" topol.top topology file potentially remove group for RES")
    parser.add_argument("RES", help=" residue name for which lines should be removed")
    parser.add_argument("OUTGRO", help=" paht to new modified gro file")
    parser.add_argument("OUTTOP", help=" path to new topolo.top file")
    args = parser.parse_args()

    # read .gro file linewise
    gro_lines = None
    with open(args.INGRO) as f:
        gro_lines = f.readlines()
    #print(gro_lines[1].strip())
    #print(len(gro_lines))
    #print('')
    # remove lines with RES
    gro_lines = [l for l in gro_lines if args.RES not in l]

    #print(gro_lines[1].strip())
    #print(len(gro_lines))
    #print('')
    
    # adjust file count in .gro file
    gro_lines[1] = str(len(gro_lines)-3) + "\n"

    #print(gro_lines[1].strip())
    #print(len(gro_lines))
    #print('')
    # write out new gro file
    with open(args.OUTGRO, "w") as f_out:
        for l in gro_lines:
            f_out.write(l)

    # read topology for potential adjusting
    # looking for a line with res in the 
    # molecules section:
    #
    # [ molecules ]
    # ; Compound       #mols
    # HOH              14881
    #
    # WARNGIN: going to assume molecules section at end of file
    # meaning to be deleted line somweher between [ molecules ] and EOF

    top_lines = []
    with open(args.INTOP) as f_top:
        top_lines = f_top.readlines()


    # WARNING: no error checking here
    ind_mol = [idx for idx, l in enumerate(top_lines) if "[ molecules ]" in l][0]
    #print(ind_mol)
    N_top_lines = len(top_lines)

    # assume RES only appears once...
    for i in range(ind_mol, N_top_lines):
        try:
            if args.RES in top_lines[i]:
                del(top_lines[i])
                break
        except:
            print(i)
            print(N_top_lines)
            exit()

    # write out new topology
    with open(args.OUTTOP, "w") as f_top_out:
        for l in top_lines:
            f_top_out.write(l)
