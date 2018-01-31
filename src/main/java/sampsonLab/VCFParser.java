package sampsonLab;

import com.google.common.collect.SetMultimap;
import com.google.errorprone.annotations.Var;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.CommonInfo;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.File;
import java.util.*;


/**
 * Created by dfermin on 11/30/16.
 */
public class VCFParser {

    static File inputVCF = null, inputVCF_tabix = null;
    //static ArrayList<VariantContext> vcList = null;

    // Constructor
    VCFParser(File vcf, File tbi) {
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
    // Function returns true if any of the values in globalFunctions.allowedSites are not set to 1
    public Boolean checkForUnreportedAllowedSites() {

        for(String s : globalFunctions.allowedSites.keySet()) {
            int status = globalFunctions.allowedSites.get(s);
            if(status == 0) return true;
        }
        return false;
    }

    /*****************************************************************************************************************/
    // Function Parses the VCF file keeping only the 'allowedSites' variants
    // No filter is applied to these cases.
    public void parseByAllowedSites() {

        // This command only works if you have the tabix tbi file for the input VCF
        VCFFileReader vcfr = new VCFFileReader(inputVCF, inputVCF_tabix);

        for(String AS : globalFunctions.allowedSites.keySet()) {
            if(globalFunctions.allowedSites.get(AS) == 1) continue;

            String chrom = AS.split(":")[0];
            int ASpos = Integer.valueOf( AS.split(":")[1] );

            CloseableIterator<VariantContext> it = vcfr.query(chrom, ASpos, ASpos);

            while(it.hasNext()) {
                VariantContext vc = it.next();
                String chr = vc.getContig();
                int pos = vc.getEnd();
                double svmProb = Double.NaN;

                CommonInfo CI = vc.getCommonInfo();
                String geneInfo = CI.getAttributeAsString("DBSNP147_GENEINFO", "#NULL");


                if( vc.getID().equalsIgnoreCase("SVM_PROBABILITY") )
                    svmProb = Double.parseDouble((String) vc.getAttribute("SVM_PROBABILITY"));

                String ref = vc.getReference().getDisplayString(); // get the reference Allele NT
                String alt = vc.getAltAlleleWithHighestAlleleCount().getDisplayString(); // get the alternative Allele NT
                String dbsnp_id = vc.getID();
                VariantInfo VI = new VariantInfo(chr, pos, dbsnp_id, ref, alt, "AS"); // AS = Allowed Site
                VI.setSvmProb(svmProb);
                VI.add(vc);
                VI.allowedSite = true;

                if(CI.hasAttribute("EFF")) {
                    VI.EFF = new EFF_Features(CI.getAttributeAsString("EFF", "#NULL"));
                }

                // Iterate over the features in 'globals::featureSet'
                // Record all of these features for the current 'VI' object
                for(String s : globalFunctions.featureSet) {

                    if(s.equalsIgnoreCase("EFF")) {
                        for(String curTS : VI.EFF.eff.keySet()) {
                            VI.fetchFeature(s, vc, curTS);
                        }
                    }
                    else if( !VI.fetchFeature(s, vc, null) ) {
                        System.err.println("\nERROR with variant: " + VI.getID() + ":\nTrying to get " + s + "\n");
                        System.exit(1);
                    }
                }

                if(VI.hasCandidateSubjects("AS", null)) {   // AS = Allowed Sites

                    if(null != VI.EFF) {
                        for(String curTS : VI.EFF.eff.keySet()) {
                            VI.printSummaryString(geneInfo, curTS);
                        }
                    }
                    else VI.printSummaryString(geneInfo, null);
                }
            }

        }
    }

    /*****************************************************************************************************************/
    // Parse the VCF file, keeping only the variant calls that overlap with the exons of our genes of interest and meet our filtering
    public void parseByExon(SetMultimap<String, Transcript> geneMap, String filter, String filterType) {

        // This command only works if you have the tabix tbi file for the input VCF
        VCFFileReader vcfr = new VCFFileReader(inputVCF, inputVCF_tabix);

        ArrayList<String> geneOrder = sortGenes(geneMap);

        int FLANK = 2; // Allow for 2 basepairs of deviation around a variant site
        for(String geneId: geneOrder) { // Iterate over the genes in the given geneMap

            ArrayList<VariantInfo> viList = new ArrayList<VariantInfo>(); // prep for this iteration

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

                        VI.passesFilter(filter, curTS.getTranscriptID()); // Determine if this variant passes our filter

                        if(VI.passedFilter) {
                            VI.addTranscript(curTS.getTranscriptID(), curTS.getGeneName());
                            viList.add(VI);
                        }


//                        if(VI.passedFilter) {
//                            // The current variant passed JEXL filtering, now only report it if there is at least one sample harboring this variant.
//                            if(VI.hasCandidateSubjects(filterType)) {
//                                VI.printSummaryString(geneId, curTS.getTranscriptID());
//                            }
//                        }
                    }
                } // end loop over exons
            } // end loop over transcripts

            reportVariants(viList, filterType);
            viList.clear();

        } // end loop over genes
    }


    /*****************************************************************************************************************/
    // Parse the VCF file, keeping only the variant calls that overlap with our transcripts of interest and meet our filtering
    public void parseByTranscript(SetMultimap<String, Transcript> geneMap, String filter, String filterType) {

        // This command only works if you have the tabix tbi file for the input VCF
        VCFFileReader vcfr = new VCFFileReader(inputVCF, inputVCF_tabix);

        ArrayList<String> geneOrder = sortGenes(geneMap);

        int FLANK = 20;
        for(String geneId: geneOrder) { // Iterate over the genes in the given geneMap

            ArrayList<VariantInfo> viList = new ArrayList<VariantInfo>(); // prep for this iteration

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

                    if(VI.passedFilter) {
                        VI.addTranscript(curTS.getTranscriptID(), curTS.getGeneName());
                        viList.add(VI);
                    }

//                    if (VI.passedFilter) {
//                        // The current variant passed filtering, now only report it if there is at least one sample harboring this variant.
//                        if (VI.hasCandidateSubjects(filterType)) {
//                            VI.printSummaryString(geneId, curTS.getTranscriptID());
//                        }
//                    }
                }
            } // end loop over transcripts

            reportVariants(viList, filterType);
            viList.clear();

        } // end loop over geneIds
    }



    /*****************************************************************************************************************/
    // Function identifies which subjects for given list of variants should be reported
    void reportVariants(ArrayList<VariantInfo> VL, String filterType) {



        if(filterType == "REC") {
            // Determine how many recessive genes each subject has.
            HashMap<String, Integer> numRecVariants = new HashMap<>();
            for(VariantInfo curVariant : VL) {
                for(String k : curVariant.genotypeMap.keySet()) numRecVariants.put(k, 0);
            }

            for(VariantInfo curVariant : VL) {
                for(String k : curVariant.genotypeMap.keySet()) {
                    Genotype_Feature g = curVariant.genotypeMap.get(k);
                    if(g.genotypeInt > 0) {
                        int n = numRecVariants.get(k) + 1;
                        numRecVariants.put(k, n);
                    }
                }
            }

            for(VariantInfo curVariant : VL) {
                boolean status = curVariant.hasCandidateSubjects("REC", numRecVariants);

                if(status) {
                    for(String curTS : curVariant.affectedTranscripts.keySet()) {
                        String gene = curVariant.affectedTranscripts.get(curTS);
                        curVariant.printSummaryString(gene, curTS);
                    }
                }
            }
        }

        // Compound heterozygosity doesn't apply to the dominant filter so 'numRecVariants is 'null'
        if(filterType == "DOM") {
            for(VariantInfo curVariant : VL) {
                boolean status = curVariant.hasCandidateSubjects("DOM", null);

                if(status) {
                    for(String curTS : curVariant.affectedTranscripts.keySet()) {
                        String gene = curVariant.affectedTranscripts.get(curTS);
                        curVariant.printSummaryString(gene, curTS);
                    }
                }
            }
        }

        if(filterType == "AS") {
            for(VariantInfo curVariant : VL) {
                boolean status = curVariant.hasCandidateSubjects("AS", null);

                if(status) {
                    for(String curTS : curVariant.affectedTranscripts.keySet()) {
                        String gene = curVariant.affectedTranscripts.get(curTS);
                        curVariant.printSummaryString(gene, curTS);
                    }
                }
            }
        }

    }





    /*****************************************************************************************************************/
    // Function sorts the gene in the given map so they occur in genomic order
    private ArrayList<String> sortGenes(SetMultimap<String, Transcript> geneMap) {
        ArrayList<String> ret = new ArrayList<>();

        List<CoordClass> listCoord = new ArrayList<CoordClass>();
        Set<String> geneSet = new HashSet<>();

        for(String gid : geneMap.keys()) {

            if( geneSet.contains(gid) ) continue; // we've recorded this gene, move on


            for(Transcript curTS : geneMap.get(gid)) {

                String obsChrom = curTS.getTrimmedChrom(); // removes the 'chr' part if present from chromosome name
                int geneStart = curTS.getGeneStart();

                // In order for this to work, we need to render all chromosomes as integers
                // This code takes care of chromosome X, Y, and MT
                int chrom = 0;
                if(obsChrom.matches("^\\d+$")) chrom = Integer.valueOf(obsChrom);
                else if(obsChrom.equalsIgnoreCase("X")) chrom = 23;
                else if(obsChrom.equalsIgnoreCase("Y")) chrom = 24;
                else if(obsChrom.matches("[MT]")) chrom = 25;
                else continue;



                listCoord.add(new CoordClass(chrom, gid, geneStart));
                geneSet.add(gid);
            }
        }


        // Sort the listCoord array by Chromosome then position
        Collections.sort(listCoord, new CoordClassComparator(
                new CoordClassChromComparator(),
                new CoordClassPositionComparator()
        ));

        // Record the genes in genomic order
        for(CoordClass C : listCoord) {
            ret.add(C.getGeneID());
        }

        return ret;
    }


}
