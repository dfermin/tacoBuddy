package sampsonLab;

//import com.google.common.base.Joiner;
import com.google.common.base.Joiner;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

import javax.script.ScriptException;
import java.util.*;
import java.util.regex.Pattern;

/**
 * Created by dfermin on 1/5/17.
 */


public class VariantInfo {
    public String chr;
    public int pos;
    public String REF;
    public String ALT;
    public String snp_id;
    public String dbsnp_id;
    public boolean allowedSite; // true means that this site bypasses all filters
    public boolean passedFilter;
    public String modelType; // REC or DOM

    public double svmProb;
    public double sample_AF; // Allele frequency for the samples

    public HashMap<String, Genotype_Feature> genotypeMap; // k = sampleID, v = genotype_feature object for this subject
    public SortedSet<String> candPatients;
    public HashMap<String, String> affectedTranscripts; // k = transcriptID, v = parent gene of transcript

    public HashMap<String, Object> userFeatures; // k = feature name, v = object (String, double, int, etc..)

    // These object only get initialized if the user decides to use them
    public ESP_Features ESP = null;
    public EFF_Features EFF = null;



    public VariantInfo(String chr, int pos, String dbID, String ref, String alt, String mt) {
        this.chr = chr;
        this.pos = pos;
        this.REF = ref;
        this.ALT = alt;
        this.dbsnp_id = dbID;
        this.modelType = mt;
        this.svmProb = 0.0;
        this.passedFilter = false;
        this.sample_AF = -1.0;
        this.allowedSite = false;

        // construct snp_id string
        snp_id = chr + ":" + String.valueOf(pos);

        // Holds the features the user wants to filter data on
        userFeatures = new HashMap<String, Object>();

        ESP = new ESP_Features();

        genotypeMap = new HashMap<String, Genotype_Feature>();
        candPatients = new TreeSet<String>();
        affectedTranscripts = new HashMap<>();


        this.checkAllowedSites();
    }

    public String getID() { return snp_id; }
    public void setSvmProb(double d) { this.svmProb = d; }


    /*********************************************************************************************/
    // Function checks to see if this variant is among the variants the user gave in the 'allowedSites' fields of the input file.
    private void checkAllowedSites() {
        String search_str = chr + ":" + String.valueOf(pos);

        if( null != globalFunctions.allowedSites ) {
            if( globalFunctions.allowedSites.keySet().contains(search_str)) {
                this.allowedSite = true;
                globalFunctions.allowedSites.put(search_str, 1); // mark this allowedSite as "observed"
            }
        }
    }


    /*********************************************************************************************/
    // Record which transcripts this variant has an impact on
    public void addTranscript(String tsID, String gene) {
        affectedTranscripts.put(tsID, gene);
    }


    /*********************************************************************************************/
    // Function records the requested data from the attributes of the variant object
    // The features are stored as objects in the 'userFeatures' map
    public boolean fetchFeature(String feat, VariantContext vc, String transcriptID) {


        // You only need the transcriptID if the user requested 'EFF' as a filter criteria

        if(feat.startsWith("ESP_")) {

            String aa = vc.getAttributeAsString("ESP_AA_AC", ".").replaceAll("\\[", "").replaceAll("\\]", "");
            String ea = vc.getAttributeAsString("ESP_EA_AC", ".").replaceAll("\\[", "").replaceAll("\\]", "");
            String tac = vc.getAttributeAsString("ESP_TAC", ".").replaceAll("\\[", "").replaceAll("\\]", "");

            ESP.ESP_AA_AC = ESP.calcMAF(aa);
            ESP.ESP_EA_AC = ESP.calcMAF(ea);
            ESP.ESP_MAX_AA_EA = Math.max(ESP.ESP_AA_AC, ESP.ESP_EA_AC);
            ESP.ESP_MAF = ESP.calcMAF(tac);

            this.userFeatures.put("ESP_AA_AC", Double.valueOf(ESP.ESP_AA_AC));
            this.userFeatures.put("ESP_EA_AC",Double.valueOf( ESP.ESP_EA_AC));
            this.userFeatures.put("ESP_MAX_AA_EA", Double.valueOf(ESP.ESP_MAX_AA_EA));

            return true;
        }
        else if(feat.equalsIgnoreCase("POLYPHEN2_HVAR")) {
            String tmp = "";
            Object o = vc.getAttribute("dbNSFP_Polyphen2_HVAR_pred", "#NULL");
            String o_type = o.getClass().getSimpleName();

            if(o_type.equalsIgnoreCase("ArrayList")) {
                tmp = globalFunctions.arrayList2String(o, ".");
            }
            else {
                tmp = vc.getAttribute("dbNSFP_Polyphen2_HVAR_pred", "#NULL").toString().replaceAll("[\\[\\]\\s]+", "");
            }

            this.userFeatures.put(feat, tmp);
        }
        else if(feat.startsWith("GERP")) {
            String key = "dbNSFP_GERP___RS";
            if(vc.hasAttribute("dbNSFP_GERP++_RS")) key = "dbNSFP_GERP++_RS";
            String tmp = "";
            Object o = vc.getAttribute(key, "#NULL");
            String o_type = o.getClass().getSimpleName();


            if(o_type.equalsIgnoreCase("ArrayList")) {
                tmp = globalFunctions.arrayList2String(o, ".");
            }
            else {
                tmp = vc.getAttribute(key, "NaN").toString().replaceAll("[\\[\\]\\s]+", "");
            }
            String[] ary = tmp.split(",");
            this.userFeatures.put(feat, Double.valueOf(ary[0]));
        }
        else if(feat.equalsIgnoreCase("MUTATIONTASTER")) {
            String tmp = "";
            Object o = vc.getAttribute("dbNSFP_MutationTaster_pred", "#NULL");
            String o_type = o.getClass().getSimpleName();

            if(o_type.equalsIgnoreCase("ArrayList")) {
                tmp = globalFunctions.arrayList2String(o, ".");
            }
            else {
                tmp = vc.getAttribute("dbNSFP_MutationTaster_pred", "#NULL").toString().replaceAll("[\\[\\]\\s]+", "");
            }
            this.userFeatures.put(feat, tmp);
        }
        else if(feat.equalsIgnoreCase("SIFT")){
            String tmp = "";
            Object o = vc.getAttribute("dbNSFP_SIFT_pred", "#NULL");
            String o_type = o.getClass().getSimpleName();

            if(o_type.equalsIgnoreCase("ArrayList")) {
                tmp = globalFunctions.arrayList2String(o, ".");
            }
            else {
                tmp = vc.getAttributeAsString("dbNSFP_SIFT_pred", "#NULL").toString().replaceAll("[\\[\\]\\s]+", "");
            }
            this.userFeatures.put(feat,tmp);
        }
        else if(feat.equalsIgnoreCase("EFF")) {
            String tmp = vc.getAttributeAsString("EFF", "#NULL");

            if(!tmp.equalsIgnoreCase("#NULL")) {
                EFF = new EFF_Features(tmp);
                tmp = transcriptID.replaceAll("\\.\\d+$", "");

                this.userFeatures.put(feat, EFF.findTS(tmp));
            }
        }
        else if(feat.equalsIgnoreCase("EXAC_MAX_MAF")) {
            double v = get_EXAC_MAX_MAF(vc);
            if( v == -1 ) this.userFeatures.put(feat, "#NULL");
            else this.userFeatures.put(feat, dbl2str(v, 4, false));
        }
        else if(feat.equalsIgnoreCase("IS_LOF")) {

            if(null == this.EFF) {
                String tmp = vc.getAttributeAsString("EFF", "#NULL");
                if (tmp != "#NULL") EFF = new EFF_Features(tmp);
            }

            String ts = transcriptID.replaceAll("\\.\\d+$", "");
            String out = EFF.findTS(ts);
            if(out.startsWith("LOF")) {
                this.userFeatures.put(feat, true);
            }
            else this.userFeatures.put(feat, false);
        }
        else {
            // If you got here then the user's selected feature is not a known "special" case listed above
            if( vc.hasAttribute(feat) ) {
                String tmp = vc.getAttributeAsString(feat, "#NULL");
                String v = catagorizeObject(tmp);
                this.userFeatures.put(feat, v);
            }
            else { this.userFeatures.put(feat, "#NULL"); }
        }


        return true; // if you got here then you had no problems in this function
    }


    /*********************************************************************************************/
    // Function computes the maximum minor allele frequency as reported by the EXAC fields
    private double get_EXAC_MAX_MAF(VariantContext vc) {
        double ret = -1;

        String popn[] = { "AFR", "AMR", "EAS", "FIN", "NFE", "SAS", "OTH" };
        ArrayList<Double> obsAF = new ArrayList<>();
        for(String curPop : popn) {
            String ac_k = "EXAC_AC_" + curPop;
            String ac_tmp = vc.getAttributeAsString(ac_k, "#NULL");

            String an_k = "EXAC_AN_" + curPop;
            String an_tmp = vc.getAttributeAsString(an_k, "#NULL");

            // If you are missing either AC or AN you can't compute AF so skip this population
            if( (an_tmp == "#NULL") || (ac_tmp == "#NULL") ) continue;

            // We are making the assumption that *IF* there are multiple values for this
            // entry, the first one list is the smallest value.
            ac_tmp = ac_tmp.replaceAll("[\\[\\]]", "");
            an_tmp = an_tmp.replaceAll("[\\[\\]]", "");
            double AC = Double.valueOf( ac_tmp.split(",")[0] );
            double AN = Double.valueOf( an_tmp.split(",")[0] );
            double MAF = AC / AN;

            if(MAF > 0.5) {
                double d = 1 - MAF;
                MAF = d;
            }

            if( !Double.isNaN(MAF) ) obsAF.add(MAF);
        }

        if(obsAF.size() > 0) {
            ret = Collections.max(obsAF); // get the largest MAF observed
        }

        return ret;
    }



    /*********************************************************************************************/
    // Function applies the user's filter to this variant. Calling this function assumes VariantInfo::fetchFeature()
    // ran successfully.
    public void passesFilter(String filter_str, String curTS) {

        boolean retVal = false;

        // if you don't give the program a filter string, all variants will pass
        if( (null == filter_str) || (filter_str.length() == 0) ) {
            passedFilter = true;
            return;
        }


        // check to see if this current variant is among the sites the user specified as 'allowedSites'
        if(this.allowedSite) {
            passedFilter = true;
            retVal = true;
        }
        else {

            if(this.userFeatures.containsKey("EFF")) {

                // Quick check to see if this is a synonymous mutation. If it is, automatically fail it
                if(EFF.isSynonymousVariant(curTS)) {
                    passedFilter = false;
                    return;
                }

                // Check to see if this is a high-impact mutation
                // Keep any high-impact mutations that have a minor allele frequency < req_min_sample_maf in our VCF file
                if(EFF.checkEFFimpact(curTS, "HIGH") && (this.sample_AF < globalFunctions.required_min_sample_af) ) {
                    passedFilter = false;
                    return;
                }
            }

            VariantFiltering VF = new VariantFiltering(filter_str);
            try {
                if(filter_str.contains("SAMPLE_AF"))
                    userFeatures.put("SAMPLE_AF", this.sample_AF);

                retVal = VF.evalVariant(userFeatures);
            } catch (ScriptException e) {
                e.printStackTrace();
            }

            if (retVal) passedFilter = true;
        }
    }

    /*********************************************************************************************/
    // Function takes in an 'object' string and determines what type of object it really is.
    // This function is used for dealing with what is returned by vc.getAttribute()
    private String catagorizeObject(String s) {
        String ret = s; // by default, just return the original value of 's' if it's a string

        // Check if 's' is 100% decimal value
        if( java.util.regex.Pattern.matches("^\\d+\\.\\d+$", s) ) {
            ret = dbl2str(Double.parseDouble(s), 4, false);
        }

        // Check if 's' is 100% integer value
        else if( java.util.regex.Pattern.matches("[^\\.\\D]+", s) ) {
            ret = Integer.toString( Integer.parseInt(s) ); // this way we get rid of '00006' cases (for example)
        }

        // Check if 's' is an array
        else if( java.util.regex.Pattern.matches("^[\\[\\(].+", s) ) {
            String tmp = s.replaceAll("[\\[\\]]", "");
            String ary[] = tmp.trim().split(",");
            ret = catagorizeObject(ary[0]); // when we encounter an array we will take the first element only
        }

        return ret;
    }


    /**********************************************************************************************/
    public void printSummaryString (String geneId, String transcriptID) {
        int PRECISION = 3;

        if(null == geneId) geneId = ""; // to handle the 'AllowedSites' instances
        if(null == transcriptID) transcriptID = ""; // to handle the 'AllowedSites' instances

        if(this.allowedSite) {
            if( !this.modelType.contains("*") )
                this.modelType += "*";
        }

        for(String k : this.candPatients) {
            Genotype_Feature gf = genotypeMap.get(k);
            ArrayList<String> ary = new ArrayList<String>();

            ary.add(k);
            ary.add(gf.genotype_word);
            ary.add(this.modelType);
            ary.add(transcriptID);
            ary.add(geneId);
            ary.add(this.chr);
            ary.add(String.valueOf(this.pos));
            ary.add(this.dbsnp_id);
            ary.add(this.REF);
            ary.add(this.ALT);
            ary.add(gf.getReadCount()); // read depth
            ary.add( dbl2str(this.sample_AF, PRECISION, true) );

            for(String feat : globalFunctions.featureSet) {
                if(this.userFeatures.containsKey(feat)) {

                    Object o = this.userFeatures.get(feat);
                    String o_type = o.getClass().getSimpleName(); // get the object type

                    if (o_type.equalsIgnoreCase("double")) {
                        double d = (Double) o;
                        if(d == 1000) ary.add(".");
                        else ary.add(dbl2str(d, PRECISION, false));
                    }
                    else if(o_type.equalsIgnoreCase("Boolean")) {
                        ary.add(String.valueOf(o));
                    }
                    else if(o_type.equalsIgnoreCase("String")) {
                        String tmp = (String) this.userFeatures.get(feat);
                        ary.add(tmp);
                    }
                    else  if(o_type.equalsIgnoreCase("ArrayList")) {
                        // Here the feature is returned as a ArrayList Object.
                        // We need to make it non-redundant and then into a string
                        String tmp = globalFunctions.arrayList2String(o, "#NULL");
                        ary.add( tmp );

                    }
                }
                else ary.add("#NULL");
            }

            String line = Joiner.on("\t").join(ary);
            System.out.println(line.replaceAll("#NULL", "."));

            ary.clear();
        }

    }

    /*********************************************************************************************/
    // Function adds the observed genotype information for all the subjects for the current variant
    public void add(VariantContext vc) {

        ArrayList<String> sampleNames = (ArrayList<String>) vc.getSampleNamesOrderedByName();
        for(String s : sampleNames) {
            Genotype G = vc.getGenotype(s);
            Genotype_Feature gf = new Genotype_Feature(s, G);
            genotypeMap.put(s, gf);
        }

        this.calcSampleAF(); // compute the sample allele frequency now that you have recorded all the samples
    }


    /*********************************************************************************************/
    // Function returns true if at least 1 subject in the genotypeMap meets the filtering requirements for this variant
    public boolean hasCandidateSubjects(String filterType, HashMap<String, Integer> numRecVariants) {
        boolean ret = false;

        for(Genotype_Feature g : genotypeMap.values()) {

            if(g.genotype_word.equalsIgnoreCase("#NULL")) continue; // no information for this variant in this sample so skip it.

            // The filter is for dominant variants
            // Keep subjects that are HET or HOM_ALT
            if( filterType.equalsIgnoreCase("DOM") && (g.genotypeInt > 0) ) candPatients.add(g.sampleID);

            // The filter is for recessive variants
            // Keep subjects that are HOM_ALT **OR**
            // Who are HET and have at least 2 recessive variants that are HOM_ALT or HET
            if( filterType.equalsIgnoreCase("REC") ) {

                int numREC_het = numRecVariants.get(g.sampleID);

                if(g.genotypeInt == 2) candPatients.add(g.sampleID);
                else if( (g.genotypeInt == 1) && (numREC_het > 1) ) candPatients.add(g.sampleID);
            }
            //if( filterType.equalsIgnoreCase("REC") && (g.genotypeInt == 2) ) candPatients.add(g.sampleID);


            // This is an allowed site outside
            if( filterType.equalsIgnoreCase("AS") && this.allowedSite && (g.genotypeInt > 0) ) candPatients.add(g.sampleID);

        }

        if(candPatients.size() > 0) ret = true;

        return ret;
    }



    /*********************************************************************************************/
    // Function returns the allele frequency (AF) for all the subjects in the VCF file.
    // The formula is: sum_over_i(genotype[0,1,2]) / (2*N)
    // Where i = a subject, and N = number of subjects
    public void calcSampleAF() {

        double sum = 0;
        double n = 0;

        for(Genotype_Feature g : genotypeMap.values()) {
            if(g.genotype_DNA_str.equalsIgnoreCase(".")) continue; // the variant wasn't called for this person so skip this instance
            n++;
            sum += g.genotypeInt;
        }

        double N = 2.0 * n;
        sample_AF = (sum / N);
    }



    /*********************************************************************************************/
    // Manually round the given double to a specific number of digits and return a String with it's value
    private static String dbl2str(double d, int precision, boolean returnPercent) {

        if(Double.isNaN(d)) return "#NULL"; // quick and dirty fix to 'NaN' values

        if(returnPercent) d *= 100.0;

        String t = "%." + Integer.toString(precision) + "f";

        String compStr = "0.";
        for(int i = 0; i < precision; i++) compStr += "0";

        // If after truncation 't' is zero, report 0 instead
        String ret = String.format(t,d);
        if(ret.equalsIgnoreCase(compStr)) ret = "0";

        return ret;
    }
}
