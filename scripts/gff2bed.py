import sys

if len(sys.argv) < 2:
    print "Usage: %s gff_file1 [gff_file2 [gff_file3...]]" % sys.argv[0]

for f_name in sys.argv[1:]:
    print f_name
    with open(f_name) as f:
        with open(f_name+".bed", "wb") as g:
            for line in f:
                if line.startswith('#'):
                    pass
                elif line.strip() == "":
                    pass
                else:
                    fields = line.split()
                    chrm = fields[0]
                    if chrm[:3] != "chr":
                        chrm = "chr" + chrm
                    start = fields[3]
                    end = fields[4]
                    sig = fields[5]
                    pk_name = fields[7]
                    nl = "\t".join([chrm, start, end, pk_name, sig]) + "\n"
                    g.write(nl)