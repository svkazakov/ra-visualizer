
import operator, random, argparse
from pyx import *
from utils import *


# ---------------------------  Global parameters  ---------------------------
zoom_LEVEL      = 2000       # 2000 nucs. in 1 point
# zoom_LEVEL         = 5000


# --------------------------  Fixed parameters  -----------------------------
reads_Y_FROM   = 3.5       + 0.3
contigs_Y      = 2.5       + 0.3
ref_Y          = 2         + 0.3
coords_bar_Y   = 1.5       + 0.3

reads_ROWS_COUNT = 50 + 1
reads_ROWS_DY    = 0.2

draw_FROM      = 3.0
legend_FROM    = 0.6

point2x_COEF   = 20     # 20 points in 1 cm coordinates    (real 35 in 1 cm)

reads_row_TO   = [11, 30, 50]
reads_row_FROM = [ 0, 15, 35]
reads_GAP_between = [10000, 2000, 100]      # in nucs.
read_NAMES = ["long", "middle-size", "short"]


class Visualizer(object):
    def __init__(self, ref_size, contigs_FN=None, reads_FN=None, view_pos_from=0, view_size=6000000, output_FN="output.svg"):
        self.ref_SIZE = ref_size
        self.contigs_FN = contigs_FN
        self.reads_FN = reads_FN
        
        self.view_POS_FROM = view_pos_from
        self.view_SIZE = view_size
        self.output_FN = output_FN
        
        self.view_POS_TO = self.view_POS_FROM + self.view_SIZE

        self.c = None        # Canvas


    # ---------------------------  Low-level functions  -----------------------------
    def coord2x(self, pos):
        assert self.view_POS_FROM <= pos <= self.view_POS_TO
        return draw_FROM + (pos - self.view_POS_FROM) / zoom_LEVEL / point2x_COEF
    
    def normalize_coords(self, pos_from, pos_to):
        start_truncated, end_truncated = False, False
        if pos_from < self.view_POS_FROM:
            start_truncated = True
            pos_from = self.view_POS_FROM
        if pos_to > self.view_POS_TO:
            end_truncated = True
            pos_to = self.view_POS_TO
        return self.coord2x(pos_from), self.coord2x(pos_to), start_truncated, end_truncated
    
    def _draw_ref(self, pos_from, pos_to, y, type="var2", size_coef=1.0):
        x0, x1, start_truncated, end_truncated = self.normalize_coords(pos_from, pos_to)
    
        self._do_draw_ref(x0, x1, y, type, size_coef, start_truncated, end_truncated)
    
    def _do_draw_ref(self, x0, x1, y, type="var2", size_coef=1.0, start_truncated=False, end_truncated=False):
        assert x0 <= x1
    #    print("in _draw_ref():   x0 = " + str(x0) + ", x1 = " + str(x1) + ", y = " + str(y))
    #    print("start_truncated = " + str(start_truncated) + ", end_truncated = " + str(end_truncated))
        line_width = style.linewidth(style._defaultlinewidth * 16 * size_coef)
        if (type == "var1") or (x1 - x0 < 0.075*2):
            self.c.stroke(path.line(x0, y, x1, y), [line_width])         # line variant 1
            return
    
        dx = 0.075
        # dx = 0
        if x0 + dx <= x1 - dx:
            x0_mod = x0 + dx
            x1_mod = x1 - dx
        else:
            x0_mod = x1_mod = (x0 + x1) / 2
    #    print("x0_mod = " + str(x0_mod) + ", x1_mod = " + str(x1_mod))
    
        self.c.stroke(path.line(x0_mod if not start_truncated else x0, y,
                           x1_mod if not end_truncated else x1, y), [line_width])       # line variant 2
        if not start_truncated:
            graph.style._diamondsymbol(self.c, unit.topt(x0_mod), unit.topt(y), unit.topt(0.17), [deco.filled([color.rgb.black])])
        if not end_truncated:
            graph.style._diamondsymbol(self.c, unit.topt(x1_mod), unit.topt(y), unit.topt(0.17), [deco.filled([color.rgb.black])])
    
    
    def _draw_read(self, pos_from, pos_to, row, styles_to_add=[], rect=False):
        x0, x1, start_truncated, end_truncated = self.normalize_coords(pos_from, pos_to)
        y = reads_Y_FROM + reads_ROWS_DY * row
    
        if rect:
            self.c.stroke(path.rect(x0, y-0.1, x1-x0, 0.2), styles_to_add)  # variant 1
            # c.stroke(path.rect(x0, y-0.1, x1, y+0.1), [style.linewidth.Thick]+styleAdd)  # variant 1
        else:
            self.c.stroke(path.line(x0, y, x1, y), [style.linewidth.THICK] + styles_to_add)  # variant 1


    # ----------------------    Main drawing functions   ------------------------
    def draw_reference(self):
        print("Drawing reference...")
        self.c.text(legend_FROM, ref_Y-0.12, r"\textbf{Reference}", [text.size.large])
        self._draw_ref(0, self.ref_SIZE - 1, ref_Y)
    
    def draw_coords_bar(self):
        print("Drawing coordinates bar...")
        self.c.stroke(path.line(self.coord2x(self.view_POS_FROM), coords_bar_Y,
                                self.coord2x(self.view_POS_TO), coords_bar_Y))
    
        for x in range(self.view_POS_FROM, self.view_POS_TO + 1, 100000):
            size = 0.05
            if x % 50000 == 0:
                size = 0.08
            if x % 100000 == 0:
                size = 0.1
    
            self.c.stroke(path.line(self.coord2x(x), coords_bar_Y + 0.01,
                                    self.coord2x(x), coords_bar_Y - size))
            if x % 1000000 == 0:
                self.c.text(self.coord2x(x) + (0.08 if x != 0 else 0), coords_bar_Y - 0.4, compressInt(x), [text.size.footnotesize, text.halign.boxcenter])
            elif x % 50000 == 0:
                self.c.text(self.coord2x(x) - 0.2, coords_bar_Y - 0.3, compressInt(x), [text.size.tiny])
    
    
    def draw_contigs(self):
        if not self.contigs_FN:
            print("Contigs file not set!  Skipping drawing contigs...")
            return
    
        self.c.text(legend_FROM, contigs_Y-0.12, r"\textbf{Contigs}", [text.size.large])
    
        print("Loading contigs from " + self.contigs_FN + " and drawing them...")
        found = 0; inside = 0
        with open(self.contigs_FN, "r") as f:
            f.readline()    # header
            last_end = None
            for line in f:
                record = line.split("\t")
                found += 1
                
                if record[2] != 'u':    # i.e. read is mapped
                    ref_start = int(record[4])
                    ref_end   = int(record[5])
                    assert 0 <= ref_start <= ref_end < self.ref_SIZE
        
                    good = True
                    if (ref_end < self.view_POS_FROM) or (ref_start > self.view_POS_TO):
                        good = False
                    else:
                        inside += 1
                        self._draw_ref(ref_start, ref_end, contigs_Y)
                        xStyle = [style.linewidth.Thin, color.cmyk.Gray, style.linestyle.dashed]
                        if self.view_POS_FROM <= ref_start <= self.view_POS_TO:
                            self.c.stroke(path.line(self.coord2x(ref_start), contigs_Y,
                                                    self.coord2x(ref_start), 14), xStyle)
                        if self.view_POS_FROM <= ref_end <= self.view_POS_TO:
                            self.c.stroke(path.line(self.coord2x(ref_end), contigs_Y,
                                                    self.coord2x(ref_end), 14), xStyle)
        
                        if last_end is not None:
                            gap_middle = (last_end + ref_start) // 2
                            gap_size = ref_start - last_end - 1
                            gap_middle_x = self.coord2x(gap_middle)
                            text_y = coords_bar_Y-1
                            path_y = text_y+0.3
                            if gap_size in [6715, 10029, 16, 12, -55, 1122]:        # move to lower position
                                text_y -= 0.3
                                path_y = text_y+0.25
                            if gap_size in [1145, 93, 5038]:
                                path_y = text_y+0.25
                            self.c.stroke(path.line(gap_middle_x, path_y, gap_middle_x, contigs_Y-0.18),
                                            [style.linewidth.Thin, color.rgb.red, deco.earrow.small])
                            self.c.text(gap_middle_x, text_y, r"\textbf{Gap= " + int2str(gap_size, delim=",") + "}",
                                            [text.size.scriptsize, color.rgb.red, text.halign.boxcenter])
                        last_end = ref_end
    
        print("Done")
        print("Total " + str(found) + " contigs found,  " + str(inside) + " are inside view window!")
        # draw_ref(     0,  20000, contigs_Y)
        # draw_ref( 25000, 100000, contigs_Y)
        # draw_ref(120000, 140000, contigs_Y)
    

    def draw_reads(self):
        if not self.reads_FN:
            print("Reads file not set!  Skipping drawing reads...")
            return
    
        self.c.text(legend_FROM, reads_Y_FROM + reads_ROWS_COUNT*reads_ROWS_DY/2 - 0.12, r"\textbf{Reads}", [text.size.large])
        # 13,  32-33
        x_style = [style.linewidth.Thin, color.cmyk.Lavender, style.linestyle.dashed]
        y = reads_Y_FROM + reads_ROWS_DY * 13
        self.c.stroke(path.line(self.coord2x(self.view_POS_FROM), y,
                                self.coord2x(self.view_POS_TO), y), x_style)
        y = reads_Y_FROM + reads_ROWS_DY * 32.5
        self.c.stroke(path.line(self.coord2x(self.view_POS_FROM), y,
                                self.coord2x(self.view_POS_TO), y), x_style)
    
    
        print("Loading reads from " + self.reads_FN + "...")
    
        all_reads = 0; good_count = 0; too_small = 0
        reads = [[], [], []]
    
        all_records = 0; using = 0
        cnt = 0
    
        with open(self.reads_FN, "r") as f:
            f.readline()    # header
            for line in f:
                record = line.strip().split('\t')
                all_reads += 1
    
                if record[2] != 'u':    # i.e. read is mapped
                    read = Read(record, line.strip(), self.ref_SIZE)
    
                    if (read.ref_end < self.view_POS_FROM) or (read.ref_start > self.view_POS_TO):
                        pass    # outside the view window
                    else:
                        # a good one
                        good_count += 1
                        if read.mapping_len >= 40000:
                            gr = 0
                        elif read.mapping_len >= 20000:
                            gr = 1
                        else:
                            gr = 2
    
                        # gr = 0      # temporary!!
                        all_records += 1
                        if read.mapping_len >= 200:
                            using += 1
                            reads[gr].append(read)
                            if read.mapping_len < 1000:
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
        print("Middle-size reads     >= 20 kb :     " + int2str(len(reads[1]),fl=5) + "     %.0f%%" % (len(reads[1])*100.0/good_count))
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
                    if read.ref_start - left >= avail_from[y]:
                        # xStyle = [color.rgb.green] if read.useful else []
    
                        # if read.ref_len < 500:
                        #     draw_read(read.ref_start-left, read.ref_end+right, y, xStyle + [style.linewidth.Thin, style.linestyle.dashed], rect=True)
                        # else:
                        self._draw_read(read.ref_start - left, read.ref_end + right, y, x_style, rect=True)
    
                        self._draw_read(read.ref_start, read.ref_end, y, x_style)
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
                        print("Warning!  Can't print " + read_NAMES[gr] + " read with length = " + int2str(read.mapping_len))
    
        print("Done!")
        print("not printed short reads = " + int2str(gr2_not_printed))

    
    def run(self):
        print("Starting...")
    
        text.set(text.LatexRunner)
        text.preamble(r"\usepackage{times}")
    
        self.c = canvas.canvas()
        self.c.stroke(path.rect(0, 0, self.coord2x(self.view_POS_TO) + 1.0, 15), [color.rgb.white])
    
        # drawing...
        self.draw_reference()
        self.draw_coords_bar()
    
        self.draw_contigs()
        self.draw_reads()
    
    
        print("Saving to svg...")
        self.c.writeSVGfile(self.output_FN)
    
        print("Done!")
        print("Resulting picture is saved to " + self.output_FN)


class Read():
    def __init__(self, r, line, ref_size):
        # r = line.split("\t")
        self.line = line

        read_name   = r[0]
        read_len    = int(r[1])
        map_strain  = r[2]
        assert map_strain in ["+", "-"]
        ref_name    = r[3]
        ref_start   = int(r[4])
        ref_end     = int(r[5])
        mapping_len = int(r[6])
        interesting = (r[7]  == "YES")
        useful      = (r[8] == "YES")
        # clipped_head= int(r[11])      # ToDo:  add to generator
        # clipped_tail= int(r[12])

        assert 0 <= ref_start <= ref_end < ref_size
        assert mapping_len >= 0

        self.read_name = read_name
        self.read_len = read_len
        self.map_strain = map_strain
        self.ref_name = ref_name
        self.ref_start= ref_start
        self.ref_end  = ref_end
        self.mapping_len  = mapping_len
        self.interesting = interesting
        self.useful = useful
        self.clipped_head = 0
        self.clipped_tail = 0



def main():
    parser = argparse.ArgumentParser(description='ra-visualizer -- Python program for visualizing genome reads and/or '
                                                 'contigs alignments to reference.')

    parser.add_argument("--ref-size", help="Reference size (default: 4'641'652 (E.coli))", type=int, default=4641652)
    parser.add_argument("-c", "--contigs", help="Contigs .stats file (generated by gen_stats_table.py)")
    parser.add_argument("-r", "--reads", help="Reads .stats file (generated by gen_stats_table.py)")

    parser.add_argument("-f", "--from-pos", help="View position: from coordinate (default: 0)", type=int, default=0)
    parser.add_argument("-s", "--view-size", help="View window size (default: 6'000'000)", type=int, default=6000000)
    parser.add_argument("-o", "--output", help="Output .svg filename (default: output.svg)", type=str, default="output.svg")

    args = parser.parse_args()
    Visualizer(args.ref_size, args.contigs, args.reads, args.from_pos, args.view_size, args.output)\
        .run()


if __name__ == '__main__':
    main()
