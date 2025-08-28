import argparse


if __name__ == "__main__":


    parser = argparse.ArgumentParser(description='Rename HOH to SOL')
    #parser.add_argument('-f','--foo', help='Description for foo argument', required=True)
    #parser.add_argument('-b','--bar', help='Description for bar argument', required=True)
    parser.add_argument("INTOP", help=" input topol.top ")
    parser.add_argument("OUTTOP", help=" output topol.top ")
    args = parser.parse_args()


    top_lines = None
    with open(args.INTOP) as f:
        top_lines = f.readlines()    


    def f_rename(l):
        
        if "HOH" in l:
            return l.replace("HOH", "SOL")
        else:
            return l 

    top_lines = [f_rename(l) for l in top_lines]

    with open(args.OUTTOP, "w") as f_out_top:
        for l in top_lines:
            f_out_top.write(l)


