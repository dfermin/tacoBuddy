package sampsonLab;

import com.google.common.collect.SetMultimap;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.File;
import java.util.ArrayList;


/**
 * Created by dfermin on 11/30/16.
 */
public class VCFParser {

    static File inputVCF = null, inputVCF_tabix = null;
    static ArrayList<VariantContext> vcList = null;

    // Constructor
    VCFParser(File vcf, File tbi) {
        inputVCF = vcf;
        inputVCF_tabix = tbi;

        if(!inputVCF_tabix.exists()) {
            System.err.print("\nERROR! Unable to find TABIX file " + inputVCF_tabix.getName() + "\n\n");
            System.exit(0);
        }
    }


    // Parse the VCF file, keeping only the variant calls that overlap with our genes of interest and meet our filtering
    public void parse(SetMultimap<String, Transcript> geneMap, String filter, String filterType) {

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
                        double svmProb = Double.NaN;
                        if( vc.getID().equalsIgnoreCase("SVM_PROBABILITY") )
                            svmProb = Double.parseDouble((String) vc.getAttribute("SVM_PROBABILITY"));

                        String ref = vc.getReference().getDisplayString(); // get the reference Allele NT
                        String alt = vc.getAltAlleleWithHighestAlleleCount().getDisplayString(); // get the alternative Allele NT


                        VariantInfo VI = new VariantInfo(chr, pos, ref, alt);
                        VI.setSvmProb(svmProb);

                        for(String s : globalFunctions.featureSet) {
                            if( !VI.fetchFeature(s, vc) ) {
                                System.err.println("\nERROR with variant: " + VI.getID() + ":\nTrying to get " + s + "\n");
                                System.exit(1);
                            }
                        }

                        // Apply the jexl filter to this variant. If it passes the
                        if( VI.filter(filter) ) { // keep this variant because it met all of our filtering criteria.
                            VI.add(vc);
                            VI.calcSampleAF();

                            if(VI.hasCandidateSubjects(filterType)) {
                                ArrayList<String> patientGenotype = VI.getPatientData();
                                System.out.print(VI.returnSummaryString(geneId, patientGenotype, filterType));
                            }
                        }

                        //System.err.print(chr + ":" + pos + ref + ">" + alt + "\n"); // progress indicator
                    }
                }
            }
        }
    }
}
