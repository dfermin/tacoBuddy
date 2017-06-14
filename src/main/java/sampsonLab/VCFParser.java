package sampsonLab;

import com.google.common.collect.SetMultimap;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.File;
import java.sql.SQLException;
import java.util.*;


/**
 * Created by dfermin on 11/30/16.
 */
public class VCFParser {

    static File inputVCF = null, inputVCF_tabix = null;
    static ArrayList<VariantContext> vcList = null;

    // Constructor
    VCFParser(File vcf, File tbi) throws SQLException {
        inputVCF = vcf;
        inputVCF_tabix = tbi;

        if(!inputVCF_tabix.exists()) {
            System.err.print("\nERROR! Unable to find TABIX file " + inputVCF_tabix.getName() + "\n\n");
            System.exit(0);
        }
    }

    /*****************************************************************************************************************/
    // Function prints out all of the unique INFO fields the user can query in the VCF file
    public void printINFOfields() {
        TreeMap<String, String> infoMap = new TreeMap<String, String>();
        VCFFileReader vcfr = new VCFFileReader(inputVCF, inputVCF_tabix);

        CloseableIterator<VariantContext> it = vcfr.iterator();
        while(it.hasNext()) {
            VariantContext vc = it.next();
            for(String k : vc.getAttributes().keySet()) {
                Object o = vc.getAttribute(k);
                infoMap.put(k, o.getClass().getSimpleName());
            }
        }
        vcfr.close();

        System.out.print("\n#INFO field\tdata type\n");
        for(String k : infoMap.keySet()) {
            System.out.println(k + "\t" + infoMap.get(k));
        }
    }


    /*****************************************************************************************************************/
    // Function returns all of the unique INFO fields the user can query in the VCF file
    public Set<String> getINFOfields() {
        HashSet<String> info = new HashSet<String>();
        VCFFileReader vcfr = new VCFFileReader(inputVCF, inputVCF_tabix);

        CloseableIterator<VariantContext> it = vcfr.iterator();
        while(it.hasNext()) {
            VariantContext vc = it.next();
            for(String k : vc.getAttributes().keySet()) {
                Object o = vc.getAttribute(k);
                info.add(k);
            }
        }
        vcfr.close();
        return info;
    }



    /*****************************************************************************************************************/
    // Parse the VCF file, keeping only the variant calls that overlap with the exons of our genes of interest and meet our filtering
    public void parseByExon(SetMultimap<String, Transcript> geneMap, String filter, String filterType) throws SQLException {

        vcList = new ArrayList<VariantContext>();

        // This command only works if you have the tabix tbi file for the input VCF
        VCFFileReader vcfr = new VCFFileReader(inputVCF, inputVCF_tabix);

        int FLANK = 2; // Allow for 2 basepairs of deviation around a variant site
        for(String geneId: geneMap.keySet()) { // Iterate over the genes in the given geneMap

            for(Transcript curTS : geneMap.get(geneId)) { // Iterate over the transcripts for this gene

                for(Exon curE : curTS.getAllExons()) { // Iterate over the exons for this transcript

                    // Get out all variant calls from the VCF file that overlap with any of the
                    // exons for this transcript
                    CloseableIterator<VariantContext> it = vcfr.query(
                            curTS.getTrimmedChrom(),
                            (curE.exonStart - FLANK),
                            (curE.exonEnd + FLANK));

                    while(it.hasNext()) {
                        VariantContext vc = it.next();
                        String chr = vc.getContig();
                        int pos = vc.getEnd();
                        double svmProb = Double.NaN;

                        if( vc.getID().equalsIgnoreCase("SVM_PROBABILITY") )
                            svmProb = Double.parseDouble((String) vc.getAttribute("SVM_PROBABILITY"));

                        String ref = vc.getReference().getDisplayString(); // get the reference Allele NT
                        String alt = vc.getAltAlleleWithHighestAlleleCount().getDisplayString(); // get the alternative Allele NT
                        String dbsnp_id = vc.getID();
                        VariantInfo VI = new VariantInfo(chr, pos, dbsnp_id, ref, alt, filterType);
                        VI.setSvmProb(svmProb);
                        VI.add(vc);


                        // Iterate over the features in 'globals::featureSet'
                        // Record all of these features for the current 'VI' object
                        for(String s : globalFunctions.featureSet) {
                            if( !VI.fetchFeature(s, vc, curTS.getTranscriptID()) ) {
                                System.err.println("\nERROR with variant: " + VI.getID() + ":\nTrying to get " + s + "\n");
                                System.exit(1);
                            }
                        }

                        // Load each variant into the database object
//                        tacoBuddy.DB.loadVariant(VI, geneId);

                        VI.passesFilter(filter, curTS.getTranscriptID()); // Determine if this variant passes our filter

                        if(VI.passedFilter) {
                            // The current variant passed JEXL filtering, now only report it if there is at least one sample harboring this variant.
                            if(VI.hasCandidateSubjects(filterType)) {
                                VI.printSummaryString(geneId, curTS.getTranscriptID());
                            }
                        }
                    }
                }
            }
        }
    }


    /*****************************************************************************************************************/
    // Parse the VCF file, keeping only the variant calls that overlap with our transcripts of interest and meet our filtering
    public void parseByTranscript(SetMultimap<String, Transcript> geneMap, String filter, String filterType) {

        vcList = new ArrayList<VariantContext>();

        // This command only works if you have the tabix tbi file for the input VCF
        VCFFileReader vcfr = new VCFFileReader(inputVCF, inputVCF_tabix);

        int ctr = 1;
        int FLANK = 20;
        for(String geneId: geneMap.keySet()) { // Iterate over the genes in the given geneMap

            for(Transcript curTS : geneMap.get(geneId)) { // Iterate over the transcripts for this gene

                // Get out all variant calls from the VCF file that overlap with this transcript
                CloseableIterator<VariantContext> it = vcfr.query(
                        curTS.getTrimmedChrom(),
                        (curTS.get_ts_Start() - FLANK),
                        (curTS.get_ts_End() + FLANK));

                while(it.hasNext()) {
                    VariantContext vc = it.next();
                    String chr = vc.getContig();
                    int pos = vc.getStart();
                    double svmProb = Double.NaN;

                    if (vc.getID().equalsIgnoreCase("SVM_PROBABILITY"))
                        svmProb = Double.parseDouble((String) vc.getAttribute("SVM_PROBABILITY"));

                    String ref = vc.getReference().getDisplayString(); // get the reference Allele NT
                    String alt = vc.getAltAlleleWithHighestAlleleCount().getDisplayString(); // get the alternative Allele NT
                    String dbsnp_id = vc.getID();

                    VariantInfo VI = new VariantInfo(chr, pos, dbsnp_id, ref, alt, filterType);
                    VI.setSvmProb(svmProb);
                    VI.add(vc);

                    // Iterate over the features in 'featureSet'
                    // Record all of these features for the current 'VI' object
                    for (String s : globalFunctions.featureSet) {
                        if (!VI.fetchFeature(s, vc, curTS.getTranscriptID())) {
                            System.err.println("\nERROR with variant: " + VI.getID() + ":\nTrying to get " + s + "\n");
                            System.exit(1);
                        }
                    }

                    VI.passesFilter(filter, curTS.getTranscriptID()); // Determine if this variant passes our filter

                    if (VI.passedFilter) {
                        // The current variant passed JEXL filtering, now only report it if there is at least one sample harboring this variant.
                        if (VI.hasCandidateSubjects(filterType)) {
                            VI.printSummaryString(geneId, curTS.getTranscriptID());
                        }
                    }
                }
            }
        }
    }
}
