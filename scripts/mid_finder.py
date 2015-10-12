import sys

if len(sys.argv) < 2:
    print "Usage: %s bed_file1 [bed_file2 [bed_file3...]]" % sys.argv[0]

for f_name in sys.argv[1:]:
    print f_name
    with open(f_name) as f:
        with open(f_name+".mid", "wb") as g:
            for line in f:
                if not line.startswith("chr"):
                    pass
                else:
                    fields = line.split()
                    start = eval(fields[1])
                    end = eval(fields[2])
                    c = fields[0]
                    pk = fields[3]
                    val = fields[4]
                    mid = (end-start)//2 + start
                    mid1 = mid + 1
                    nl = "\t".join([c, str(mid), str(mid1), pk, val]) + "\n"
                    g.write(nl)