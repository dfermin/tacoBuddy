package sampsonLab;

import com.google.common.base.Joiner;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import javafx.beans.binding.BooleanExpression;
import org.apache.commons.jexl3.*;

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
    public int allowedSite; // -1 = not allowed, 0 = DOM, 1 = REC
    public boolean passedFilter;
    public String modelType; // REC or DOM

    public double svmProb;
    public double sample_MAF; // minor allele frequency for the samples

    public HashMap<String, Genotype_Feature> genotypeMap; // k = sampleID, v = genotype_feature object for this subject
    public SortedSet<String> candPatients;

    public HashMap<String, Object> userFeatures; // k = feature name, v = object (String, double, int, etc..)

    // These object only get initialized if the user decides to use them
    public ESP_Features ESP = null;
    public EFF_Features EFF = null;
    public dbNSFP_Features dbNSFP = null;



    public VariantInfo(String chr, int pos, String dbID, String ref, String alt, String mt) {
        this.chr = chr;
        this.pos = pos;
        this.REF = ref;
        this.ALT = alt;
        this.dbsnp_id = dbID;
        this.modelType = mt;
        this.svmProb = 0.0;
        this.passedFilter = false;
        this.sample_MAF = -1.0;
        this.allowedSite = -1;

        // construct snp_id string
        snp_id = chr + ":" + String.valueOf(pos);

        userFeatures = new HashMap<String, Object>();

        ESP = new ESP_Features();
        dbNSFP = new dbNSFP_Features();

        genotypeMap = new HashMap<String, Genotype_Feature>();
        candPatients = new TreeSet<String>();

        this.checkAllowedSites();
    }


    public String getID() { return snp_id; }
    public void setSvmProb(double d) { this.svmProb = d; }

    /*********************************************************************************************/
    // Function checks to see if this variant is among the variants the user gave in the 'allowedSites' fields of the input file.
    private void checkAllowedSites() {
        String search_str = chr + ":" + String.valueOf(pos);

        // user gave no allowedSites input.
        if( null == globalFunctions.allowedSitesMap ) {
            this.allowedSite = -1;
            return;
        }

        for(String k : globalFunctions.allowedSitesMap.keySet()) {
            if(search_str.equalsIgnoreCase(k)) {
                if(globalFunctions.allowedSitesMap.get(k).equalsIgnoreCase("DOM")) this.allowedSite = 0;
                else if(globalFunctions.allowedSitesMap.get(k).equalsIgnoreCase("REC")) this.allowedSite = 1;
                break;
            }
        }
    }



    /*********************************************************************************************/
    // Function records the requested data from the attributes of the variant object
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

            userFeatures.put("ESP_AA_AC", ESP.ESP_AA_AC);
            userFeatures.put("ESP_EA_AC", ESP.ESP_EA_AC);
            userFeatures.put("ESP_MAX_AA_EA",ESP.ESP_MAX_AA_EA);

            return true;
        }

        if(feat.equalsIgnoreCase("POLYPHEN2_HVAR")) {
            String tmp = vc.getAttribute("dbNSFP_Polyphen2_HVAR_pred", "#NULL").toString().replaceAll("[\\[\\]\\s]+", "");
            userFeatures.put(feat, tmp);
        }

        if(feat.equalsIgnoreCase("GERP")) {
            String key = "dbNFSP_GERP___RS";
            if(vc.hasAttribute("dbNSFP_GERP++_RS")) key = "dbNSFP_GERP++_RS";

            String tmp = vc.getAttribute(key, "NaN").toString().replaceAll("[\\[\\]\\s]+", "");
            String[] ary = tmp.split(",");
            userFeatures.put(feat, Double.valueOf(ary[0]));
        }

        if(feat.equalsIgnoreCase("MUTATIONTASTER")) {
            userFeatures.put(feat, vc.getAttribute("dbNSFP_MutationTaster_pred", "#NULL"));
        }

        if(feat.equalsIgnoreCase("SIFT")){
            String tmp = vc.getAttributeAsString("dbNSFP_SIFT_pred", "#NULL").toString().replaceAll("[\\[\\]\\s]+", "");
            userFeatures.put(feat,tmp);
        }

        if(feat.equalsIgnoreCase("EFF")) {
            String tmp = vc.getAttributeAsString("EFF", "#NULL");
            if(!tmp.equalsIgnoreCase("#NULL")) {
                EFF = new EFF_Features(tmp);
                tmp = transcriptID.replaceAll("\\.\\d+$", "");
                userFeatures.put(feat, EFF.findTS(tmp));
            }
        }

        if(feat.equalsIgnoreCase("IS_LOF")) {
            if( vc.hasAttribute("LOF") ) {
                String tmp = vc.getAttributeAsString("LOF", "#NULL");
                if (tmp.equalsIgnoreCase("#NULL")) userFeatures.put(feat, false);
                else userFeatures.put(feat, Boolean.valueOf(tmp));
            }
            else userFeatures.put(feat, false);
        }


        return true; // if you got here then you had no problems in this function
    }


    /*********************************************************************************************/
    // Function applies the user's filter to this variant. Calling this function assumes VariantInfo::fetchFeature()
    // ran successfully.
    public boolean passesFilter(String jexl_filter_str, String curTS) {

        boolean retVal = false;

        // check to see if this current variant is among the sites the user specified as 'allowedSites'
        if(this.allowedSite > -1) {
            passedFilter = true;
            retVal = true;
        }
        else {

            if(this.userFeatures.containsKey("EFF")) {

                // Quick check to see if this is a synonymous mutation. If it is, automatically fail it
                if(EFF.isSynonymousVariant(curTS)) return  false;

                // Check to see if this is a high-impact mutation
                // Keep any high-impact mutations that have a minor allele frequency < req_min_sample_maf in our VCF file
                if(EFF.checkEFFimpact(curTS, "HIGH") && (this.sample_MAF < globalFunctions.required_min_sample_maf) ) return true;
            }

//            System.out.println(jexl_filter_str);

            JexlEngine jexl = new JexlBuilder().cache(512).strict(true).silent(false).create(); // Create a jexl engine
            JexlExpression expr = jexl.createExpression(jexl_filter_str); // define the expression you want to test/use

            // Create a MapContext object and populate it with the variables that are in VariantInfo objects
            JexlContext jc = new MapContext();

            jc.set("SAMPLE_MAF", sample_MAF);

            for(String k : this.userFeatures.keySet()) {
                Object o = this.userFeatures.get(k);
                String dataType = o.getClass().getSimpleName();

                // The user may have selected to report a feature they are not filtering the data on.
                // This if-statement prevents an error with JEXL in these cases
                if( !jexl_filter_str.toUpperCase().contains(k.toUpperCase()) ) continue;

                if(dataType.equalsIgnoreCase("Double")) {
                    jc.set(k, o);
                }

                if(dataType.equalsIgnoreCase("String")) {
                    String tmp = (String) o;
                    jc.set(k, formatStringForJexl(tmp.replaceAll("#NULL", "#")));
                }

                if(dataType.equalsIgnoreCase("Boolean")) {
                    jc.set(k, o);
                }
            }

            retVal = (Boolean) expr.evaluate(jc);

            if (retVal) passedFilter = true;
        }

        return retVal;
    }


    /*********************************************************************************************/
    // Function formats the given string so that Jexl will view it as an array.
    // This function assumes the data passed to it is comma separated
    private String formatStringForJexl(String inputStr) {
        String ret = "";
        HashSet<String> S = new HashSet<String>();
        for(String s : inputStr.split(",")) {
            String tmp = "'" + s.trim() + "'";
            S.add(tmp);
        }
        ret = "[" + Joiner.on(",").join(S) + "]";
        return ret;
    }


    /*********************************************************************************************/
    // Function returns a non-redundant representation of the given string
    private String makeStringNR(String s) {
        String ret = formatStringForJexl(s);
        return(ret.replaceAll("[\\[\\]\\s']+", ""));
    }


    /**********************************************************************************************/
    public void printSummaryString (String geneId, String transcriptID) {
        int PRECISION = 3;

        if(this.allowedSite != -1) this.modelType += "*";

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
            ary.add( dbl2str(this.sample_MAF, PRECISION, true) );

            for(String feat : globalFunctions.featureSet) {
                if(this.userFeatures.containsKey(feat)) {

                    Object o = this.userFeatures.get(feat);

                    if (o.getClass().getSimpleName().equalsIgnoreCase("double")) {
                        double d = (Double) o;
                        if(d == 1000) ary.add(".");
                        else ary.add(dbl2str(d, PRECISION, false));
                    }
                    else if(o.getClass().getSimpleName().equalsIgnoreCase("Boolean")) {
                        ary.add(String.valueOf(o));
                    }
                    else {
                        String tmp = (String) this.userFeatures.get(feat);
                        ary.add( this.makeStringNR(tmp) );
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
    // Function adds the observed genotype information for all the subject for the current variant
    public void add(VariantContext vc) {

        ArrayList<String> sampleNames = (ArrayList<String>) vc.getSampleNamesOrderedByName();
        for(String s : sampleNames) {
            Genotype G = vc.getGenotype(s);
            Genotype_Feature gf = new Genotype_Feature(s, G);
            genotypeMap.put(s, gf);
        }

        this.calcSampleAF();
    }


    /*********************************************************************************************/
    // Function returns true if at least 1 subject in the genotypeMap meets the filtering requirements for this variant
    public boolean hasCandidateSubjects(String filterType) {
        boolean ret = false;

        for(Genotype_Feature g : genotypeMap.values()) {

            if(g.genotype_word.equalsIgnoreCase("#NULL")) continue; // no information for this variant in this sample so skip it.

            // This is a variant that is not in 'allowedSites' and the current subject is homozygous dominant so skip it..
            if( this.allowedSite < 0 && g.genotype_word.equalsIgnoreCase("HOM") ) continue;

            // This is a variant in the Dominant 'allowedSites' category but the current subject is homozygous recessive so skip it..
            if( this.allowedSite == 0 && g.genotype_word.equalsIgnoreCase("HOM_ALT") ) continue;


            // Only report 'allowedSites' from the Recessive category if the subject is HOM_ALT for this variant
            if( this.allowedSite == 1 && !g.genotype_word.equalsIgnoreCase("HOM_ALT") ) continue;

            candPatients.add(g.sampleID);
        }

        if(candPatients.size() > 0) ret = true;

        return ret;
    }

    /*********************************************************************************************/
    // Function returns the minor allele frequency (MAF) for all the subjects in the VCF file.
    // The formula is: 1 - sum_over_i(genotype[0,1,2]) / (2*N)
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
        double tmp = (sum / N);
        if(Double.isNaN(tmp)) sample_MAF = 0;
        else sample_MAF = tmp;
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
