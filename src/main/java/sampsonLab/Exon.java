package sampsonLab;

/**
 * Created by dfermin on 12/2/16.
 */



public class Exon implements Comparable<Exon> {
    String exonID;
    String geneName;
    String chromosome;
    int exonNumber, exonStart, exonEnd;

    public Exon(String eid, String gene_id, String chr, int s, int e, int i) {
        exonID = eid;
        geneName = gene_id;
        chromosome = chr;
        exonStart = s;
        exonEnd = e;
        exonNumber = i;
    }

    @Override
    public boolean equals(Object E) {
        boolean ret = false;
        if(E == null || E.getClass() != this.getClass() ) ret = false;
        if(E == this) ret = true;

        Exon Other = (Exon) E;
        if( this.exonID.equalsIgnoreCase(Other.exonID) ) ret = true;

        return ret;
    }

    public int compareTo(Exon e) {
        if(this.exonStart <= e.exonStart) return this.exonStart;
        else return e.exonStart;
    }

}
