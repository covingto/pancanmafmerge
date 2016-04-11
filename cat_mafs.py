
import os, os.path, sys

def main(args):
    fpaths = args.FPATHS
    counter = 0
    with open(args.OUTPUT, 'w') as fo:
        with open(fpaths[0], 'r') as fi:
            for line in fi.readlines():
                counter += 1
                fo.write(line)
        for fnum, fpath in enumerate(fpaths[1:]):
            with open(fpath, 'r') as fi:
                # skip til after the Hugo_Symbol line
                writing = False
                for line in fi.readlines():
                    counter += 1
                    if writing:
                        fo.write(line)
                    else:
                        if 'Hugo' in line:
                            writing = True
            print "Merged %s of %s files, Wrote %s lines total" % (str(fnum + 1), str(len(fpaths)), str(counter))

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('OUTPUT', type = str, help = 'output file')
    parser.add_argument('FPATHS', nargs = '+', type = str, help = 'input files')

    args = parser.parse_args()

    main(args)
