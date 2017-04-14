package sampsonLab;

import java.util.HashSet;

/**
 * Created by dfermin on 11/30/16.
 */
public class Transcript {
    private String chromosome;
    private int geneStart;
    private int geneEnd;
    private int startPos;
    private int endPos;
    private int tsLen;
    private char strand;
    private double tims;

    // These are features from the last column of the GFF3 file.
    // At the time this code was written, were using the field names used by gencode v25, 2016-07-18
    private String geneID = null;
    private String geneType = null;
    private String geneName = null;
    private String transcriptID = null;

    private HashSet<Exon> exonsSet = null;

    public Transcript(String ch,
                      int geneS, int geneE,
                      char c,
                      String gene_id, String gene_type, String gene_name,
                      int tsS, int tsE,
                      String tsID) {
        chromosome = ch;
        geneStart = geneS;
        geneEnd = geneE;
        strand = c;
        geneID = gene_id;
        geneType = gene_type;
        geneName = gene_name;
        startPos = tsS;
        endPos = tsE;
        transcriptID = tsID;

        tims = 10000;

        exonsSet = new HashSet<Exon>();
    }


    public String getGeneName() { return geneName; }
    public String getTranscriptID() { return transcriptID; }
    public String getChrom() { return chromosome; }
    public int get_ts_Start() { return startPos; }
    public int get_ts_End() { return endPos; }
    public int get_gene_Start() { return geneStart; }
    public int get_gene_End() { return geneEnd; }
    public int getNumExons() { return exonsSet.size(); }
    public int getTsLen() { return tsLen; }
    public double getTims() { return tims; }
    public void addExon(sampsonLab.Exon e) { exonsSet.add(e); }

    public HashSet<sampsonLab.Exon> getAllExons() { return exonsSet; }

    public String getTrimmedChrom() {
        String ret = chromosome.replace("chr", "");
        return(ret);
    }

    public void setTIMS(double TIMS) {
        this.tims = TIMS;
    }

    public void calcCDSlength() {
        tsLen = 0;
        for(Exon e : exonsSet) {
            int diff = 0;
            if(e.exonEnd < e.exonStart) diff = (e.exonStart - e.exonEnd + 1);
            else diff = (e.exonEnd - e.exonStart + 1);

            tsLen += diff;
        }
    }


    @Override
    public boolean equals(Object T) {
        boolean ret = false;
        if(T == null || T.getClass() != this.getClass()) ret = false;
        if(T == this) ret = true;

        Transcript Other = (Transcript) T;
        if( this.transcriptID.equalsIgnoreCase(Other.getTranscriptID()) ) ret = true;

        return ret;
    }


}
