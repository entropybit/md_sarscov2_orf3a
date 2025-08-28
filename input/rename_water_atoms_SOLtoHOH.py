import argparse


if __name__ == "__main__":


    parser = argparse.ArgumentParser(description='Renames water atoms from OW, HW1, HW2 to O, H1, H2')
    #parser.add_argument('-f','--foo', help='Description for foo argument', required=True)
    #parser.add_argument('-b','--bar', help='Description for bar argument', required=True)
    parser.add_argument("INGRO", help=" input GROFILE.gro")
    parser.add_argument("OUTGRO", help=" paht to new modified gro file")
    args = parser.parse_args()


    gro_lines = None
    with open(args.INGRO) as f:
        gro_lines = f.readlines()

    #gro_lines = [l for l in gro_lines if args.RES not in l]


    def f_rename(l):
    
        if "SOL" in l:
            l = l.replace("SOL", "HOH")
            if "OW" in l:
                return l.replace("OW", " O")
            if "HW1" in l:
                return l.replace("HW1", " H1")
            if "HW2" in l:
                return l.replace("HW2", " H2")

        return l


    gro_lines = [f_rename(l) for l in gro_lines]

    with open(args.OUTGRO, "w") as f_out:
        for l in gro_lines:
            f_out.write(l)


