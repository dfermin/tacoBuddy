package sampsonLab;

import java.util.HashMap;

/**
 * Created by dfermin on 3/16/17.
 */
public class Gene {
    String geneID;
    String geneType;
    String geneName;
    String geneGroup; // DOM or REC

    String chrom;
    int gStart, gEnd;

    HashMap<String, Transcript> memberTS;



}


