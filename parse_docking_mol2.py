#!/bin/python
import numpy
import optparse

def main(file):
    file=open(file)
    count=0
    separate=False
    top=False
    pattern=None
    for line in file.readlines():
        if 'MOLECULE' in line:
            separate=True
            prefix=line
            continue
        if separate==True: 
            if len(line.split())>0:
                print pattern
                pattern=line.split()[0]
                top=True
                new=open('%s.mol2' % pattern, 'w')
                new.write(prefix)
                new.write(line)
                separate=False
                continue
        elif top==True:
            new.write(line)
                
def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-i', '--input', dest='file',
                      help='file to parse')
    (options, args) = parser.parse_args()
    return (options, args)

if __name__ == "__main__":
    (options, args) = parse_commandline()
    main(file=options.file)
