package sampsonLab;

import com.google.common.collect.SetMultimap;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.File;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Set;


/**
 * Created by dfermin on 11/30/16.
 */
public class VCFParser {

    static File inputVCF = null, inputVCF_tabix = null;
    static ArrayList<VariantContext> vcList = null;
    static Connection conn = null;

    VCFParser(File vcf, File tbi) {
        inputVCF = vcf;
        inputVCF_tabix = tbi;
    }

    public void parse(SetMultimap<String, Transcript> geneMap) throws SQLException {

        vcList = new ArrayList<VariantContext>();
        String q = null;
        Statement stmt = null;

        String DBfileName = "jdbc:hsqldb:file:" + (new File("VCFDB")).getAbsolutePath() + ";shutdown=true"; // just down the DB when last connection is closed

        try {
            Class.forName("org.hsqldb.jdbc.JDBCDriver");
        } catch (ClassNotFoundException e) {
            e.printStackTrace();
        }
        try {

            conn = DriverManager.getConnection(DBfileName, "SA", "");
            q = "CREATE TABLE rawd (" +
                    "gene_id VARCHAR(50), " +
                    "trans_id VARCHAR(50), " +
                    "exon_id VARCHAR(50), " +
                    "exonNum INTEGER, " +
                    "variant VARCHAR(100), " +
                    "dbSNPid VARCHAR(20), " +
                    "svmProb DOUBLE, " +
                    "sampleId VARCHAR(50), " +
                    "GT VARCHAR(10), " +
                    "DP INTEGER, " +
                    "GQ INTEGER )";

            stmt = conn.createStatement();
            stmt.executeUpdate(q);

        } catch (SQLException e) {
            e.printStackTrace();
        }

        // This command only works if you have the tabix tbi file for the input VCF
        VCFFileReader vcfr = new VCFFileReader(inputVCF, inputVCF_tabix);

        // Iterate over the genes in the given geneMap
        for(String geneId: geneMap.keySet()) {

            Set<Transcript> allTS = geneMap.get(geneId);

            for(Transcript curTS : allTS) {

                for(Exon curE : curTS.getAllExons()) {

                    // Get out all variant calls from the VCF file that overlap with any of the
                    // exons for this transcript
                    CloseableIterator<VariantContext> it = vcfr.query(
                            curTS.getTrimmedChrom(),
                            curE.exonStart,
                            curE.exonEnd);

//                    while(it.hasNext()) {
//                        VariantContext vc = it.next();
//                        vcList.add( vc );
//                    }

                    //debugging
                    while(it.hasNext()) {

                        VariantContext vc = it.next();
                        String chr = vc.getContig();
                        int pos = vc.getEnd();
                        double svmProb = Double.parseDouble((String) vc.getAttribute("SVM_PROBABILITY"));
                        String ref = vc.getReference().getDisplayString(); // get the reference Allele nucleotide
                        String alt = vc.getAltAlleleWithHighestAlleleCount().getDisplayString();
                        String dbSNP_id = vc.getID();

                        Genotype G = vc.getGenotype(0);


                        q = "INSERT INTO rawd (" +
                            " gene_id, " +
                            " trans_id, " +
                            " exon_id, " +
                            " exonNum, " +
                            " variant, " +
                            " dbSNPid, " +
                            " svmProb, " +
                            " sampleId, " +
                            " GT, DP, GQ) VALUES ( ";

                        q +=  "'" + curTS.getGeneName() + "', " +
                                "'" + curTS.getTranscriptID() + "', " +
                                "'" + curE.exonID + "', " +
                                    + curE.exonNumber + ", " +
                                "'" + chr + ":" + pos + ref + ">" + alt + "', " +
                                "'" + dbSNP_id + "', " +
                                    + svmProb + ", " +
                                "'" + G.getSampleName() + "', " +
                                "'" + G.getGenotypeString() + "', " +
                                    + G.getDP() + ", " +
                                    + G.getGQ() + ")";

                        stmt.executeUpdate(q);
                        conn.commit();
                    }
                }
            }
        }
        vcfr.close();

        stmt.close();
        conn.close();

        System.err.print("Parsed " + inputVCF.getName() + ", " + vcList.size() + " variant calls kept.\n");
    }


    public void db_variants() {

        Iterator<VariantContext> it = vcList.iterator();
        while(it.hasNext()) {
            VariantContext vc = (VariantContext) it.next();
            Genotype G = vc.getGenotype((172-10));

            System.out.print(
                    "chr" + vc.getContig() + ":" + vc.getEnd() + vc.getReference().getDisplayString() + ">" + vc.getAltAlleleWithHighestAlleleCount().getDisplayString() + "\t" +
                    G.getSampleName() + "\t" + G.getGenotypeString() + "\n");


            break;
        }
    }
}
