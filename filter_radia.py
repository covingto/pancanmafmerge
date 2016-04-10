
import os, os.path, sys

def main(args):
    with open(args.INPUT, 'r') as fi, open(args.OUTPUT, 'w') as fo:
        for line in fi.readlines():
            if line[0] != '#' and 'SS=1' in line and 'PASS' in line:
                continue
            fo.write(line)

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('INPUT', type = str, help = 'input file path')
    parser.add_argument('OUTPUT', type = str, help = 'output file path')

    args = parser.parse_args()

    main(args)

