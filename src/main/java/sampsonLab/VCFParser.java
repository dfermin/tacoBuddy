package sampsonLab;

import com.google.common.collect.SetMultimap;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import java.io.File;
import java.sql.*;
import java.util.ArrayList;
import java.util.Collection;
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
        PreparedStatement prep = null;

        String DBfileName = "jdbc:hsqldb:file:" + (new File("VCFDB")).getAbsolutePath() + ";shutdown=true"; // just down the DB when last connection is closed

        try {
            Class.forName("org.hsqldb.jdbc.JDBCDriver");
        } catch (ClassNotFoundException e) {
            e.printStackTrace();
        }
        try {

            conn = DriverManager.getConnection(DBfileName, "SA", "");
            stmt = conn.createStatement();
            stmt.executeUpdate("DROP TABLE IF EXISTS rawd");
            stmt.execute("SET FILES LOG FALSE");
            stmt.execute("SET FILES NIO SIZE 2048");

            q = "CREATE CACHED TABLE rawd (" +
                    "gene_id VARCHAR(50), " +
                    "trans_id VARCHAR(50), " +
                    "exon_id VARCHAR(50), " +
                    "exonNum INTEGER, " +
                    "variant VARCHAR(100), " +
                    "dbSNPid VARCHAR(20), " +
                    "svmProb DOUBLE, " +
                    "sampleId VARCHAR(50), " +
                    "GT VARCHAR(50), " +
                    "DP INTEGER, " +
                    "GQ INTEGER )";
            stmt.executeUpdate(q);

        } catch (SQLException e) {
            e.printStackTrace();
        }

        q = "INSERT INTO rawd (gene_id, trans_id, exon_id, exonNum, variant, dbSNPid, svmProb, sampleId, GT, DP, GQ) " +
                "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?) ";
        prep = conn.prepareStatement(q);

        // This command only works if you have the tabix tbi file for the input VCF
        VCFFileReader vcfr = new VCFFileReader(inputVCF, inputVCF_tabix);


        VCFHeader hdr = vcfr.getFileHeader();
        Collection<VCFInfoHeaderLine> H = hdr.getInfoHeaderLines();
        for(VCFInfoHeaderLine h : H) {
            System.out.println(h.getID() + "\t" + h.getDescription());
        }
        System.exit(0);



        int ctr = 1;
        for(String geneId: geneMap.keySet()) { // Iterate over the genes in the given geneMap

            Set<Transcript> allTS = geneMap.get(geneId);

            for(Transcript curTS : allTS) {

                for(Exon curE : curTS.getAllExons()) {

                    // Get out all variant calls from the VCF file that overlap with any of the
                    // exons for this transcript
                    CloseableIterator<VariantContext> it = vcfr.query(
                            curTS.getTrimmedChrom(),
                            curE.exonStart,
                            curE.exonEnd);

                    while(it.hasNext()) {
                        VariantContext vc = it.next();
                        String chr = vc.getContig();
                        int pos = vc.getEnd();
                        double svmProb = Double.parseDouble((String) vc.getAttribute("SVM_PROBABILITY"));
                        String ref = vc.getReference().getDisplayString(); // get the reference Allele nucleotide
                        String alt = vc.getAltAlleleWithHighestAlleleCount().getDisplayString();
                        String dbSNP_id = vc.getID();

                        int N = vc.getNSamples();
                        for(int i = 0; i < vc.getNSamples(); i++) { // iterate over the samples
                            Genotype G = vc.getGenotype(i);


                            String variantStr = chr + ":" + pos + ref + ">" + alt;

                            prep.setString(1, curTS.getGeneName() );
                            prep.setString(2, curTS.getTranscriptID() );
                            prep.setString(3, curE.exonID );
                            prep.setInt(4, curE.exonNumber );
                            prep.setString(5, variantStr );
                            prep.setString(6, dbSNP_id );
                            prep.setDouble(7, svmProb );
                            prep.setString(8, G.getSampleName() );
                            prep.setString(9, G.getGenotypeString() );
                            prep.setInt(10, G.getDP() );
                            prep.setInt( 11, G.getGQ() );

                            prep.addBatch();

                            if( (ctr%10000) == 0 ) {
                                conn.setAutoCommit(false);
                                prep.executeBatch();
                                conn.setAutoCommit(true);
                                prep.clearBatch();
                            }
                            ctr++;
                        }
                    }
                }
            }
        }
        vcfr.close();
        conn.commit();

        stmt.execute("SET FILES LOG TRUE");

        stmt.close();
        conn.close();
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
