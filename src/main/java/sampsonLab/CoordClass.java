package sampsonLab;

import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

/**
 * Created by dfermin on 1/29/18.
 */
public class CoordClass implements Comparable<CoordClass> {

    private int chrom;
    private String geneID;
    private int startPos;

    // Used to sort Coordinates
    public int compareTo(CoordClass a) {
        return this.startPos - a.startPos;
    }

    // Default constructor
    public CoordClass(int c_, String g_, int s_) {
        this.chrom = c_;
        this.geneID = g_;
        this.startPos = s_;
    }

    public int getChrom() { return chrom; }
    public String getGeneID() { return geneID; }
    public int getStartPos() { return startPos; }

    @Override
    public String toString() {
        String c = "";
        if(chrom == 23) c = "X";
        if(chrom == 24) c = "Y";
        if(chrom == 25) c = "M";
        else c = String.valueOf(chrom);
        return String.format("%s\t%s:%d", geneID, c, startPos);
    }
}


class CoordClassComparator implements Comparator<CoordClass> {

    private List<Comparator<CoordClass>> listComparators;

    @SafeVarargs
    public CoordClassComparator(Comparator<CoordClass>... comparators) {
        this.listComparators = Arrays.asList(comparators);
    }

    @Override
    public int compare(CoordClass a, CoordClass b) {
        for(Comparator<CoordClass> comparator : listComparators) {
            int result = comparator.compare(a, b);
            if(result != 0) {
                return result;
            }
        }
        return 0;
    }
}

class CoordClassChromComparator implements Comparator<CoordClass> {
    @Override
    public int compare(CoordClass a, CoordClass b) {
        return a.getChrom() - b.getChrom();
    }
}


class CoordClassPositionComparator implements Comparator<CoordClass> {
    @Override
    public int compare(CoordClass a, CoordClass b) {
        return a.getStartPos() - b.getStartPos();
    }
}
