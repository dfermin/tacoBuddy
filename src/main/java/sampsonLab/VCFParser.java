package sampsonLab;

import com.google.common.collect.SetMultimap;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;


/**
 * Created by dfermin on 11/30/16.
 */
public class VCFParser {

    static File inputVCF = null, inputVCF_tabix = null;
    static ArrayList<VariantContext> vcList = null;

    VCFParser(File vcf, File tbi) {
        inputVCF = vcf;
        inputVCF_tabix = tbi;
    }


    // Parse the VCF file, keeping only the variant calls that overlap with our genes of interest and meet our filtering
    public void parse(SetMultimap<String, Transcript> geneMap, String filter) {

        vcList = new ArrayList<VariantContext>();

        // This command only works if you have the tabix tbi file for the input VCF
        VCFFileReader vcfr = new VCFFileReader(inputVCF, inputVCF_tabix);


        int ctr = 1;
        for(String geneId: geneMap.keySet()) { // Iterate over the genes in the given geneMap

            for(Transcript curTS : geneMap.get(geneId)) { // Iterate over the transcripts for this gene

                for(Exon curE : curTS.getAllExons()) { // Iterate over the exons for this transcript

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
                        String ref = vc.getReference().getDisplayString(); // get the reference Allele NT
                        String alt = vc.getAltAlleleWithHighestAlleleCount().getDisplayString(); // get the alternative Allele NT

                        VariantInfo VI = new VariantInfo(chr, pos, ref, alt);

                        for(String s : globalFunctions.featureSet) {
                            if( !VI.fetchFeature(s, vc) ) {
                                System.err.println("\nERROR with variant: " + VI.getID() + ":\nTrying to get " + s + "\n");
                                System.exit(1);
                            }
                        }

                        VI.filter(filter);

                        int j = 1;
                    }
                }
            }
        }


//        int ctr = 1;
//        for(String geneId: geneMap.keySet()) { // Iterate over the genes in the given geneMap
//
//            Set<Transcript> allTS = geneMap.get(geneId);
//
//            for(Transcript curTS : allTS) {
//
//                for(Exon curE : curTS.getAllExons()) {
//
//                    // Get out all variant calls from the VCF file that overlap with any of the
//                    // exons for this transcript
//                    CloseableIterator<VariantContext> it = vcfr.query(
//                            curTS.getTrimmedChrom(),
//                            curE.exonStart,
//                            curE.exonEnd);
//
//                    while(it.hasNext()) {
//                        VariantContext vc = it.next();
//                        String chr = vc.getContig();
//                        int pos = vc.getEnd();
//                        double svmProb = Double.parseDouble((String) vc.getAttribute("SVM_PROBABILITY"));
//                        String ref = vc.getReference().getDisplayString(); // get the reference Allele NT
//                        String alt = vc.getAltAlleleWithHighestAlleleCount().getDisplayString(); // get the alternative Allele NT
//                        String dbSNP_id = vc.getID();
//
//
//
//                        int N = vc.getNSamples();
//                        for(int i = 0; i < vc.getNSamples(); i++) { // iterate over the samples
//                            Genotype G = vc.getGenotype(i);
//
//
//                            String variantStr = chr + ":" + pos + ref + ">" + alt;
//
//                            prep.setString(1, curTS.getGeneName() );
//                            prep.setString(2, curTS.getTranscriptID() );
//                            prep.setString(3, curE.exonID );
//                            prep.setInt(4, curE.exonNumber );
//                            prep.setString(5, variantStr );
//                            prep.setString(6, dbSNP_id );
//                            prep.setDouble(7, svmProb );
//                            prep.setString(8, G.getSampleName() );
//                            prep.setString(9, G.getGenotypeString() );
//                            prep.setInt(10, G.getDP() );
//                            prep.setInt( 11, G.getGQ() );
//
//                            prep.addBatch();
//
//                            if( (ctr%10000) == 0 ) {
//                                conn.setAutoCommit(false);
//                                prep.executeBatch();
//                                conn.setAutoCommit(true);
//                                prep.clearBatch();
//                            }
//                            ctr++;
//                        }
//                    }
//                }
//            }
//        }
//        vcfr.close();
//        conn.commit();
//
//        stmt.execute("SET FILES LOG TRUE");
//
//        stmt.close();
//        conn.close();
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
