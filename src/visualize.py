
import operator, random, argparse
from pyx import *
from utils import *




# ----------------------------  Global parameters  ----------------------------


view_POS_FROM      =         0
view_SIZE          =   6000000
# ref_SIZE           = 4641652
# ref_SIZE           = 133797422
ref_SIZE           = None

view_POS_TO        = view_POS_FROM + view_SIZE

zoom_LEVEL         = 2000        #  1000 nucs. in 1 point
# zoom_LEVEL         = 5000


# dataset             = "Human_chr10"
# contigs_FN          = "../../nanopore-data-server/{dataset}/edges.stats"    .format(dataset=dataset)
# reads_FN            = "../../nanopore-data-server/{dataset}/dataset.stats3"  .format(dataset=dataset)
contigs_FN = None
reads_FN = None


# --------------------------  Fixed parameters  -----------------------------


reads_Y_FROM   = 3.5        +0.3
contigs_Y      = 2.5        +0.3
ref_Y          = 2          +0.3
coords_bar_Y   = 1.5        +0.3


reads_ROWS_COUNT = 50  +1
reads_ROWS_DY    = 0.2


draw_FROM    = 3.0
legend_FROM  = 0.6

point2x_COEF = 20           #  20 points in 1 cm coordinates    (real 35 in 1 cm)



c = None        # Canvas

# ---------------------------  Low-level functions  -----------------------------


def coord2x(pos):
    assert view_POS_FROM <= pos <= view_POS_TO
    return draw_FROM + (pos - view_POS_FROM) / zoom_LEVEL / point2x_COEF

def normalizeCoords(posFrom, posTo):
    start_truncated, end_truncated = False, False
    if posFrom < view_POS_FROM:
        start_truncated = True
        posFrom = view_POS_FROM
    if posTo > view_POS_TO:
        end_truncated = True
        posTo = view_POS_TO
    return coord2x(posFrom), coord2x(posTo), start_truncated, end_truncated


def draw_ref(posFrom, posTo, y, type="var2", size_coef=1.0):
    x0, x1, start_truncated, end_truncated = normalizeCoords(posFrom, posTo)

    _draw_ref(x0, x1, y, type, size_coef, start_truncated, end_truncated)

def _draw_ref(x0, x1, y, type="var2", size_coef=1.0, start_truncated=False, end_truncated=False):
    assert x0 <= x1
#    print("in _draw_ref():   x0 = " + str(x0) + ", x1 = " + str(x1) + ", y = " + str(y))
#    print("start_truncated = " + str(start_truncated) + ", end_truncated = " + str(end_truncated))
    line_width = style.linewidth(style._defaultlinewidth * 16 * size_coef)
    if (type == "var1") or (x1 - x0 < 0.075*2):
        c.stroke(path.line(x0, y, x1, y), [line_width])         # variant 1
        return

    dx = 0.075
    # dx = 0
    if x0 + dx <= x1 - dx:
        x0_mod = x0 + dx
        x1_mod = x1 - dx
    else:
        x0_mod = x1_mod = (x0 + x1) / 2
#    print("x0_mod = " + str(x0_mod) + ", x1_mod = " + str(x1_mod))

    c.stroke(path.line(x0_mod if not start_truncated else x0, y,
                       x1_mod if not end_truncated   else x1, y), [line_width])       # variant 2
    if not start_truncated:
        graph.style._diamondsymbol(c, unit.topt(x0_mod), unit.topt(y), unit.topt(0.17), [deco.filled([color.rgb.black])])
    if not end_truncated:
        graph.style._diamondsymbol(c, unit.topt(x1_mod), unit.topt(y), unit.topt(0.17), [deco.filled([color.rgb.black])])


def draw_read(posFrom, posTo, row, styles_to_add=[], rect=False):
    x0, x1, start_truncated, end_truncated = normalizeCoords(posFrom, posTo)
    y = reads_Y_FROM + reads_ROWS_DY * row

    if rect:
        c.stroke(path.rect(x0, y-0.1, x1-x0, 0.2), styles_to_add)  # variant 1
        # c.stroke(path.rect(x0, y-0.1, x1, y+0.1), [style.linewidth.Thick]+styleAdd)  # variant 1
    else:
        c.stroke(path.line(x0, y, x1, y), [style.linewidth.THICK] + styles_to_add)  # variant 1


# ----------------------    main drawing functions   ------------------------

def draw_reference():
    print("Drawing reference...")
    c.text(legend_FROM, ref_Y-0.12, r"\textbf{Reference}", [text.size.large])
    draw_ref(0, ref_SIZE - 1, ref_Y)

def draw_coords_bar():
    print("Drawing coordinates bar...")
    c.stroke(path.line(coord2x(view_POS_FROM), coords_bar_Y, coord2x(view_POS_TO), coords_bar_Y))

    for x in range(view_POS_FROM, view_POS_TO+1, 100000):
        size = 0.05
        if x % 50000 == 0:
            size = 0.08
        if x % 100000 == 0:
            size = 0.1

        c.stroke(path.line(coord2x(x), coords_bar_Y+0.01, coord2x(x), coords_bar_Y - size))
        if x % 1000000 == 0:
            c.text(coord2x(x)+(0.08 if x != 0 else 0), coords_bar_Y-0.4, compressInt(x), [text.size.footnotesize, text.halign.boxcenter])
        elif x % 50000 == 0:
            c.text(coord2x(x)-0.2, coords_bar_Y-0.3, compressInt(x), [text.size.tiny])


def draw_contigs():
    if not contigs_FN:
        print("Contigs file not set!  Skipping drawing contigs...")
        return

    c.text(legend_FROM, contigs_Y-0.12, r"\textbf{Contigs}", [text.size.large])

    print("Loading contigs from " + contigs_FN + " and drawing them...")
    found = 0; inside = 0
    with open(contigs_FN, "r") as f:
        f.readline()    # header
        last_end = None
        for line in f:
            record = line.split("\t")
            found += 1

            ref_start = int(record[4])
            ref_end   = int(record[5])
            assert 0 <= ref_start <= ref_end < ref_SIZE

            good = True
            if (ref_end < view_POS_FROM) or (ref_start > view_POS_TO):
                good = False
            else:
                inside += 1
                draw_ref(ref_start, ref_end, contigs_Y)
                xStyle = [style.linewidth.Thin, color.cmyk.Gray, style.linestyle.dashed]
                if view_POS_FROM <= ref_start <= view_POS_TO:
                    c.stroke(path.line(coord2x(ref_start),contigs_Y, coord2x(ref_start),14), xStyle)
                if view_POS_FROM <=  ref_end  <= view_POS_TO:
                    c.stroke(path.line(coord2x(ref_end)  ,contigs_Y, coord2x(ref_end),  14), xStyle)

                if last_end is not None:
                    gap_middle = (last_end + ref_start) // 2
                    gap_size = ref_start - last_end - 1
                    gapMiddleX = coord2x(gap_middle)
                    textY = coords_bar_Y-1
                    pathY = textY+0.3
                    if gap_size in [6715, 10029, 16, 12, -55, 1122]:        # move to lower position
                        textY -= 0.3
                        pathY = textY+0.25
                    if gap_size in [1145, 93, 5038]:
                        pathY = textY+0.25
                    c.stroke(path.line(gapMiddleX, pathY, gapMiddleX, contigs_Y-0.18),
                                    [style.linewidth.Thin, color.rgb.red, deco.earrow.small])
                    c.text(gapMiddleX, textY, r"\textbf{Gap= " +int2str(gap_size, delim=",") + "}",
                                    [text.size.scriptsize, color.rgb.red, text.halign.boxcenter])
                last_end = ref_end


    print("Done")
    print("Total " + str(found) + " contigs found,  " + str(inside) + " are inside view window!")
    # draw_ref(     0,  20000, contigs_Y)
    # draw_ref( 25000, 100000, contigs_Y)
    # draw_ref(120000, 140000, contigs_Y)


reads_row_TO   = [11, 30, 50]
reads_row_FROM = [ 0, 15, 35]
reads_GAP_between = [10000, 2000, 100]        # in nucs.
read_NAMES = ["long", "middle-sized", "short"]

def draw_reads():
    if not reads_FN:
        print("Reads file not set!  Skipping drawing reads...")
        return

    c.text(legend_FROM, reads_Y_FROM + reads_ROWS_COUNT*reads_ROWS_DY/2 - 0.12, r"\textbf{Reads}", [text.size.large])
    # 13,  32-33
    x_style = [style.linewidth.Thin, color.cmyk.Lavender, style.linestyle.dashed]
    y = reads_Y_FROM + reads_ROWS_DY * 13
    c.stroke(path.line(coord2x(view_POS_FROM), y, coord2x(view_POS_TO), y), x_style)
    y = reads_Y_FROM + reads_ROWS_DY * 32.5
    c.stroke(path.line(coord2x(view_POS_FROM), y, coord2x(view_POS_TO), y), x_style)


    print("Loading reads from " + reads_FN + "...")

    all_reads = 0; good_count = 0; too_small = 0
    reads = [[], [], []]

    all_records = 0; using = 0
    cnt = 0

    with open(reads_FN, "r") as f:
        f.readline()    # header
        for line in f:
            record = line.strip().split('\t')
            all_reads += 1

            if (record[4] != "u") and (record[4] != "noInfo"):
                read = Read(record, line.strip())

                if (read.ref_end < view_POS_FROM) or (read.ref_start > view_POS_TO):
                    pass    # outside the view window
                else:
                    # a good one
                    good_count += 1
                    if read.ref_len >= 40000:
                        gr = 0
                    elif read.ref_len >= 20000:
                        gr = 1
                    else:
                        gr = 2

                    gr = 0      # temporary!!
                    all_records += 1
                    if read.ref_len >= 200:
                        using += 1
                        reads[gr].append(read)
                        if read.ref_len < 1000:
                            too_small += 1

                    # if (gr != 2) and (41800000 <= read.ref_start < 41900000):
                    #     print("read  #" + read.line)
                    #     cnt += 1

    print("Done!")
    print("Total reads found :     " + int2str(all_reads))
    print()
    # print("cnt = " + str(cnt))
    print("All records           :     " + int2str(all_records))
    print("  with match >= 200   :     " + int2str(using))
    print()


    print("Short reads            < 20 kb :     " + int2str(len(reads[2]),fl=5) + "     %.0f%%" % (len(reads[2])*100.0/good_count))
    print("Middle-sized reads    >= 20 kb :     " + int2str(len(reads[1]),fl=5) + "     %.0f%%" % (len(reads[1])*100.0/good_count))
    print("Long reads            >= 40 kb :     " + int2str(len(reads[0]),fl=5) + "     %.0f%%" % (len(reads[0])*100.0/good_count))
    print("too small  = " + str(too_small))

    print("Processing...")
    for gr in [0,1,2]:
        reads[gr].sort(key=operator.attrgetter('ref_start'))

    print("and drawing them...")
    random.seed(1234688)
    avail_from = [0] * reads_ROWS_COUNT
    gr2_not_printed = 0

    for gr in [0,1,2]:
        for read in reads[gr]:
            first_y = random.randint(reads_row_FROM[gr], reads_row_TO[gr])
            if first_y == 11:
                first_y -= 1
            y = first_y
            printed = False

            if read.map_strain == '+':
                x_style = [color.rgb.blue]
                left = read.clipped_head
                right = read.clipped_tail
            else:
                assert read.map_strain == '-'
                x_style = [color.cmyk.CarnationPink]
                left = read.clipped_tail
                right = read.clipped_head

            # if read.ref_len in [27652, 20347]:
            #     xStyle = [color.rgb.red]

            while True:
                if read.ref_start-left >= avail_from[y]:
                    # xStyle = [color.rgb.green] if read.useful else []

                    # if read.ref_len < 500:
                    #     draw_read(read.ref_start-left, read.ref_end+right, y, xStyle + [style.linewidth.Thin, style.linestyle.dashed], rect=True)
                    # else:
                    draw_read(read.ref_start-left, read.ref_end+right, y, x_style, rect=True)

                    draw_read(read.ref_start, read.ref_end, y, x_style)
                    avail_from[y] = read.ref_end + reads_GAP_between[gr] + 1 + right
                    # avail_from[y] = read.ref_end + reads_GAP_between[gr] + 1
                    printed = True
                    break
                y += 1
                if y > reads_row_TO[gr]:
                    y = reads_row_FROM[gr]
                if y == first_y:
                    break
            if not printed:
                if gr == 2:
                    gr2_not_printed += 1
                else:
                    print("Warning!  Can't print " + read_NAMES[gr] + " read with length = " + int2str(read.ref_len))

    print("Done!")
    print("not printed short reads = " + int2str(gr2_not_printed))



class Read():
    def __init__(self, r, line):
        # r = line.split("\t")
        self.line = line

        no          = int(r[0])
        map_strain  = r[4]
        assert map_strain in ["+", "-"]
        ref_name    = r[5]
        ref_start   = int(r[6])
        ref_end     = int(r[7])
        ref_len     = int(r[8])
        interesting = (r[9]  == "YES")
        useful      = (r[10] == "YES")
        clipped_head= int(r[11])
        clipped_tail= int(r[12])

        assert 0 <= ref_start <= ref_end < ref_SIZE
        assert ref_len >= 0

        self.no = no
        self.map_strain = map_strain
        self.ref_name = ref_name
        self.ref_start= ref_start
        self.ref_end  = ref_end
        self.ref_len  = ref_len
        self.interesting = interesting
        self.useful = useful
        self.clipped_head = clipped_head
        self.clipped_tail = clipped_tail




def run(args):
    print("Starting...")
    print("args = " + str(args))

    global ref_SIZE, contigs_FN, reads_FN
    ref_SIZE = args.ref_size
    contigs_FN = args.contigs
    reads_FN = args.reads


    text.set(text.LatexRunner)
    text.preamble(r"\usepackage{times}")

    global c
    c = canvas.canvas()
    c.stroke(path.rect(0, 0, coord2x(view_POS_TO) + 1.0, 15), [color.rgb.white])

    # drawing...
    draw_reference()
    draw_coords_bar()

    draw_contigs()
    draw_reads()


    #start = 38911000
    #len1 = 1655
    #draw_ref(start, start + len1 -1, contigs_Y)
    #draw_ref(start, start + len1 -1, contigs_Y+0.5, type="var1")

    # draw_ref(4, 5, contigs_Y, type="var1")
    # draw_ref(4, 5, ref_Y, type="var2")

    # c.stroke(path.line(1.09, 3, 10, 3), [style.linewidth.THICK, style.linejoin.round])
    # graph.style._diamondsymbol(c, unit.topt(1), unit.topt(3), unit.topt(0.17), [deco.filled([color.rgb.black])])
    #
    # c.stroke(path.line(7, 5, 14, 5), [style.linewidth.THICK])
    # c.text(5, 4.9, r"\textbf{Reference}", [text.size.large])


    print("Saving to svg...")
    c.writeSVGfile("test.svg")

    print("Done!")


def main():
    parser = argparse.ArgumentParser(description='ra-visualizer -- Python program for visualizing genome reads and/or '
                                                 'contigs alignments to reference.')

    # parser.add_argument("print_string", help="Prints the supplied argument.")
    parser.add_argument("-s", "--ref-size", help="Reference size (default: 4'641'652 (E.coli))", type=int, default=4641652)
    parser.add_argument("-c", "--contigs", help="Contigs .stats file (generated by sam_gen.py)")
    parser.add_argument("-r", "--reads", help="Reads .stats file (generated by sam_gen.py)")

    args = parser.parse_args()
    run(args)


if __name__ == '__main__':
    main()
