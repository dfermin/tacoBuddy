package sampsonLab;

import java.util.HashMap;

/**
 * Created by dfermin on 11/30/16.
 */
public class Transcript {
    private String chromosome;
    private int geneStart;
    private int geneEnd;
    private int startPos;
    private int endPos;
    private char strand;

    // These are features from the last column of the GFF3 file.
    // At the time this code was written, were using the field names used by gencode v25, 2016-07-18
    private String geneID = null;
    private String geneType = null;
    private String geneName = null;
    private String transcriptID = null;


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
    }


    public String getGeneName() { return geneName; }
    public String getTranscriptID() { return transcriptID; }
    public String getChrom() { return chromosome; }
    public int get_ts_Start() { return startPos; }
    public int get_ts_End() { return endPos; }
    public int get_gene_Start() { return geneStart; }
    public int get_gene_End() { return geneEnd; }

    public String getTrimmedChrom() {
        String ret = chromosome.replace("chr", "");
        return(ret);
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
