package sampsonLab;

/**
 * Created by dfermin on 11/30/16.
 */
public class Gene {
    static String chromosome;
    static int startPos;
    static int endPos;
    static char strand;
    static boolean isValid;

    // These are features from the last column of the GFF3 file.
    // At the time this code was written, were using the field names used by gencode v25, 2016-07-18
    static String geneID = null;
    static String geneType = null;
    static String geneName = null;

    // Default constructor
    public Gene(String src_line) {
        isValid = false; // set to true only if you get valid data

        String[] data = src_line.split("\t");

        if(data[2].equalsIgnoreCase("gene")) {
            chromosome = data[0];
            startPos = Integer.parseInt(data[3]);
            endPos = Integer.parseInt(data[4]);
            strand = data[6].charAt(0);

            String info[] = data[8].split(";");
            for (String s : info) {
                if (s.startsWith("ID=")) geneID = s.substring(3);
                if (s.startsWith("gene_type=")) geneType = s.substring(10);
                if (s.startsWith("gene_name=")) geneName = s.substring(10);
            }
            isValid = true;
        }
    }

//    // Copy Constructor
//    public Gene(Gene another) {
//        this.chromosome = another.chromosome;
//        this.startPos = another.startPos;
//        this.endPos = another.endPos;
//        this.strand = another.strand;
//        this.isValid = another.isValid;
//        this.geneID = another.geneID;
//        this.geneType = another.geneType;
//        this.geneName = another.geneName;
//    }

    public boolean isValidGene() { return isValid; }
    public String getGeneName() { return geneName; }
    public int getStart() { return startPos; }
    public int getEnd() { return endPos; }
}
